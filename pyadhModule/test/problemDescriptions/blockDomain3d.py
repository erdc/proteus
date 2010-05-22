#!/usr/bin/env python

import numpy

def genPoly(polyfileBase="blockDomain",
            nx=4,ny=4,nz=4,
            Lx=1.0,Ly=1.0,Lz=1.0):
    """
    create a simple block domain in 3d
    """
    dx=Lx/nx; dy = Ly/ny; dz = Lz/nz

    vertices = []
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                vertices.append((0.0+dx*i,0.0+dy*j,0.0+dz*k))
            #
        #
    #
    nVertices = len(vertices)
    assert nVertices == (nx+1)*(ny+1)*(nz+1)

    #domain is x,y,z in [left,right] x [front,back], [top,bottom]
    #
    boundaryTags = {'bottom':1,
                    'right':2,
                    'top':3,
                    'left':4,
                    'front':5,
                    'back':6,
                    'interior':0}
    
    base = 1
    #write a square for each face
    nFacetsTotal = nz*ny*(nx+1) + nz*nx*(ny+1) + nx*ny*(nz+1)
    facets = {}
    nFacets = 0
    #xfaces
    #i=0, left, i = Nx right
    for i in range(nx+1):
        for j in range(ny):
            for k in range(nz):
                #
                vbb = i + j*(nx+1) + k*(nx+1)*(ny+1);
                vbt = i + (j+1)*(ny+1) + k*(nx+1)*(ny+1);
                vtb = i + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
                vtt = i + (j+1)*(ny+1) + (k+1)*(nx+1)*(ny+1);
                tag = boundaryTags['interior']
                if i == 0: tag = boundaryTags['left']
                if i == nx: tag = boundaryTags['right']
                #segment number, start vertex, end vertex, id
                facets[nFacets] = (nFacets,vbb,vbt,vtt,vtb,tag)
                nFacets += 1
    #
    #yfaces
    #j=0, back, j=Ny, front
    for j in range(ny+1):
        for i in range(nx):
            for k in range(nz):
                vlb = i + j*(nx+1) + k*(nx+1)*(ny+1);
                vrb = i+1 + j*(nx+1) +  k*(nx+1)*(ny+1);
                vlt = i + j*(nx+1) + (k+1)*(nx+1)*(ny+1);
                vrt = i+1 + j*(nx+1) +  (k+1)*(nx+1)*(ny+1);
                tag = boundaryTags['interior']
                if j == 0: tag = boundaryTags['back']
                if j == ny: tag = boundaryTags['front']
                facets[nFacets] = (nFacets,vlb,vrb,vrt,vlt,tag)
                nFacets += 1
    #
    #zfaces
    #k=0, bottom, k=Nz, top
    for k in range(nz+1):
        for i in range(nx):
            for j in range(ny):
                vlb = i + j*(nx+1) + k*(nx+1)*(ny+1);
                vrb = i+1 + j*(nx+1) +  k*(nx+1)*(ny+1);
                vlf = i + (j+1)*(nx+1) + k*(nx+1)*(ny+1);
                vrf = i+1 + (j+1)*(nx+1) +  k*(nx+1)*(ny+1);
                tag = boundaryTags['interior']
                if k == 0: tag = boundaryTags['bottom']
                if k == nz: tag = boundaryTags['top']
                facets[nFacets] = (nFacets,vlb,vrb,vrf,vlf,tag)
                nFacets += 1
    #
    assert nFacets == nFacetsTotal

    #return a table to identify regions by a unique flag too
    regions = {}
    curRegion = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                #region number, x,y,z region id
                regions[(i,j,k)] =(curRegion + base,0.+(i+0.5)*dx,0.+(j+0.5)*dy,0.+(k+0.5)*dz,curRegion)
                curRegion += 1
            #
        #
    #
    polyfile = open(polyfileBase+'.poly','w')
    polyfile.write('#poly file for [%s,%s,%s] domain with %s x %s x %s blocks \n' % (Lx,Ly,Lz,nx,ny,nz))
    polyfile.write('%d %d %d %d \n' % (nVertices,3,0,1))
    polyfile.write('#vertices \n')
    for iv in range(len(vertices)):
        polyfile.write('%d %12.5e %12.5e %12.5e %d \n' % (iv+base,vertices[iv][0],vertices[iv][1],vertices[iv][2],1))
    #

    #write facets
    polyfile.write("#facets \n")
    polyfile.write('%d %d \n' % (nFacets,1))
    for facet in range(nFacets):
        polyfile.write('%d %d %d \n' % (1,0,facets[facet][-1]))#1 polygon 0 holes, tag
        polyfile.write('%d %d %d %d %d \n' % (4, #square, and four vertices
                                               facets[facet][1]+base,
                                               facets[facet][2]+base,
                                               facets[facet][3]+base,
                                               facets[facet][4]+base))
    polyfile.write('#holes\n 0\n')
    polyfile.write('#regions\n')
    nRegions = nx*ny*nz
    polyfile.write('%d \n' % nRegions)
    curRegion = 0
    for k in range(nz):
        for j in range(ny):
                for i in range(nx):
                    polyfile.write('%d %12.5e %12.5e %12.5e %d \n' % (regions[(i,j,k)][0],
                                                                      regions[(i,j,k)][1],
                                                                      regions[(i,j,k)][2],
                                                                      regions[(i,j,k)][3],
                                                                      regions[(i,j,k)][4]))
    #
    return (Lx,Ly,Lz),boundaryTags,regions
if __name__ == '__main__':
    Lx = 1.0; Ly = 1.0; Lz = 1.0; nx = 2; ny = 2; nz = 2

    L,bt,reg = genPoly("blockDomain3d",
                       nx=nx,ny=ny,nz=nz,
                       Lx=Lx,Ly=Ly,Lz=Lz)

    print "domain= %s nx=%s ny=%s nz= %s boundaryTags=%s regions=%s " % (L,
                                                                         nx,ny,nz,
                                                                         bt,
                                                                         reg)
