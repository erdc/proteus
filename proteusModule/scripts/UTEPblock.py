#!/usr/bin/env python
import numpy

Lx=1.0; Ly=1.0;
nx=128; ny=128;
#nx=16; ny=16;
dx=Lx/nx; dy = Ly/ny

vertices = []
for j in range(ny+1):
    for i in range(nx+1):
        vertices.append((0.0+dx*i,0.0+dy*j))
    #
#
nVertices = len(vertices)
assert nVertices == (nx+1)*(ny+1)
base = 1
polyfile = open('UTEPexample.poly','w')
polyfile.write('#poly file for UTEP example 128 x 128 \n')
polyfile.write('%d %d %d %d \n' % (nVertices,2,1,0))
polyfile.write('#vertices \n')
for iv in range(len(vertices)):
    polyfile.write('%d %12.5e %12.5e %d \n' % (iv+base,vertices[iv][0],vertices[iv][1],1))
#

#write a segment for each edge
nSegments = ny*(nx+1) + nx*(ny+1)
polyfile.write('%d %d \n' % (nSegments,1))
polyfile.write('#segments \n')
curSeg = 0
for j in range(ny+1): #horizontal edges first
    for i in range(nx):
        vl = i + j*(ny+1); vr = i+1 + j*(ny+1)
        polyfile.write('%d %d %d %d \n ' % (curSeg+base,vl+base,vr+base,1))
        curSeg += 1
    #
#
for i in range(nx+1): #vertical edges
    for j in range(ny):
        vb = i + j*(ny+1); vt = i + (j+1)*(ny+1)
        polyfile.write('%d %d %d %d \n ' % (curSeg+base,vb+base,vt+base,1))
        curSeg += 1
    #
#
polyfile.write('#holes\n 0\n')
polyfile.write('#regions\n')
nRegions = nx*ny
polyfile.write('%d \n' % nRegions)
curRegion = 0
for j in range(ny):
    for i in range(nx):
        polyfile.write('%d %12.5e %12.5e %d \n' % (curRegion+base,0.+(i+0.5)*dx,0.+(j+0.5)*dy,curRegion))
        curRegion +=1
#
