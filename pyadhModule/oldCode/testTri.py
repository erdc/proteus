#!/usr/bin/env python
#import standard Python modules
import sys,os
#module for triangle
import triangulate
import numpy as Numeric
#define some flags for how the script will run
verbose = 2   #level of output

trin = triangulate.new()
print 'trin= ',trin
pointsin =  Numeric.array([[0.0, 0.0],
                        [1.0,0.0],
                        [1.0,1.0],
                        [0.0,1.0]])
marksin  = Numeric.ones((4,),Numeric.Int)

triangulate.setPoints(trin,pointsin)

pointsout = triangulate.getPoints(trin)
pointsout2 = triangulate.getPointsCopy(trin)
marksout = triangulate.getPointMarkers(trin)
print 'pointsin= ',pointsin
print 'pointsout= ',pointsout
print 'pointsoutCopy= ',pointsout2
print 'pointmarksin= ',marksin
print 'pointmarksout= ',marksout

#create some point attributes,
npts  = pointsout.shape[0]
natts = 1
patts = Numeric.zeros((npts,natts),Numeric.Float)
for i in range(npts):
    patts[i,:] = pointsout[i,0]**2 + pointsout[i,1]**2
#end i
triangulate.setPointAttributes(trin,patts)

pattsout = triangulate.getPointAttributes(trin)
pattsoutCopy = triangulate.getPointAttributesCopy(trin)
print 'pattsin=\n',patts
print 'pattsout=\n',pattsout
print 'pattsoutCopy=\n',pattsoutCopy

#manually construct triangles
ntri = 2;
ncorn= 3;
elems0  = Numeric.zeros((ntri,ncorn),Numeric.Int)
elems0[0,:] = [0,1,2]
elems0[1,:] = [0,2,3]

triangulate.setTriangles(trin,elems0)

#shallow copy
elems1 = triangulate.getTriangles(trin)
#deep copy
elems1Copy = triangulate.getTrianglesCopy(trin)

print 'elems0=\n',elems0
print 'elems1=\n',elems1

#construct boundary segments
nsegs = 4
segsin= Numeric.zeros((nsegs,2),Numeric.Int)
segsin[0,:] = [3,0]  #left
segsin[1,:] = [1,2]  #right
segsin[2,:] = [0,1]  #bottom
segsin[3,:] = [2,3]  #top
segmarks = Numeric.zeros((nsegs,),Numeric.Int)
segmarks[:] = [0,1,2,3]
print 'setting segments=\n',segsin
print 'segmarks=\n',segmarks

triangulate.setSegments(trin,segsin)

segsout    = triangulate.getSegments(trin)
segmarksout= triangulate.getSegmentMarkers(trin)

print 'segmentsOut=\n',segsout
print 'segmarksout=\n',segmarksout


ntatts = 2
nelem  = elems1.shape[0]
tatts = Numeric.zeros((nelem,ntatts),Numeric.Float)
for i in range(nelem):
    tatts[i,:] = Numeric.ones(ntatts,)*float(i)
#end i
print 'try to set triangle attributes \n',tatts
triangulate.setTriangleAttributes(trin,tatts)

carea = Numeric.ones((nelem,),Numeric.Float)
print 'try to set triangle area constraints \n',carea
triangulate.setTriangleAreas(trin,carea)

#this should be empty before calling triangulate
neigs0 = triangulate.getNeighbors(trin)
print 'neigs= ',neigs0
#this should be empty before calling triangulate
neigs1 = triangulate.getNeighborsCopy(trin)
print 'neigsCopy= ',neigs1

vornorms= triangulate.getVoronoiNormals(trin)
print 'vornorms= ',vornorms

#get basic info about triangulation
info =triangulate.getInfo(trin)
sinfo = """triangulation:
number of points        = %d
attributes per point    = %d
number of triangles     = %d
nodes per triangle      = %d
attributes per triangle = %d
number of segments      = %d
number of holes         = %d
number of regions       = %d
number of edges         = %d
"""
print sinfo % info


#create a simple flag list
flags = "rzenqv"

triout = triangulate.new()
vorout = triangulate.new()

triangulate.applyTriangulate(flags,trin,triout,vorout)

#get basic info about triangulation output

oinfo =triangulate.getInfo(triout)
print 'output triangulation has'
print sinfo % oinfo

vinfo =triangulate.getInfo(vorout)
print 'voronio output has'
print sinfo % vinfo

print 'calling printReport for output'
triangulate.printReport(triout)


#delete trin manually to see what happens with data storage
print 'deleting trin manually'
del trin
print 'should be ok: pointsoutCopy=\n',pointsout2
print 'should be garbage pointsout=\n',pointsout
print 'should be ok: elems1Copy=\n',elems1Copy
print 'should be garbage elems1=\n',elems1
print 'should be ok: pattsoutCopy=\n',pattsoutCopy
print 'should be garbage patts=\n',pattsout

#make sure triout still has its data around
elems2 = triangulate.getTriangles(triout)
print 'triout triangles=\n',elems2
