#!/usr/bin/env  python
import numpy as numpy
#path to input file
#nx=300; ny=300;
xr=[50,250]; yr=[200,300];
nx=xr[1]-xr[0]; ny=yr[1]-yr[0];

readSliceDirectly=True#False
if readSliceDirectly == True:
    datafile2d = "E6_jun1_50x250_200x300.dat" #"./E6_jun1_slice.dat"
    f = open(datafile2d,'r')
    a = f.read()
    A = a.split()
    AA= [int(AI) for AI in A]
    #mwf debug
    #import pdb
    #pdb.set_trace()
    B=numpy.array(AA,'i')
    aSlice=numpy.reshape(B,(nx,ny))
else:
    datafile3d = "/Users/mfarthin/Public/mTex/chrisNotes/particleScale6.1/E6_jun1.bin"
    f = open(datafile3d,'rb')
    #this is 300**3 image
    a = f.read()
    a.split()
    #len(a)
    aSlice = numpy.zeros((nx,ny),'i')
    k = 150*300*2
    for i in range(xr[0],xr[1]):
        for j  in range(yr[0],yr[1]):
            #print  a[i*300+j]
            if a[k+i*300+j] =='\x00':
                aSlice[i-xr[0],j-yr[0]] = 0

            if a[k+i*300+j] =='\x01':
                aSlice[i-xr[0],j-yr[0]] = 1

            if a[k+i*300+j] =='\x02':
                aSlice[i-xr[0],j-yr[0]] = 2
        #
    #
#
dat = open("slice.dat",'w')
ncolor = [0,0,0]
for i in range(nx):
    for j  in range(ny):
        if aSlice[i,j] == 0:
            ncolor[0] += 1
        elif aSlice[i,j] == 1:
            ncolor[1] += 1
        elif aSlice[i,j] == 2:
            ncolor[2] += 1
        dat.write('%i \t' % aSlice[i,j])
    dat.write('\n')
dat.close()
#set up funcs to build R,G,and B arrays from the gray scale 0,1,2 data then use merge('RGB') to build color image

#experiment with writing to a file that triangle can mesh
L=(1.0,1.0); dx=L[0]/nx; dy=L[1]/ny;
base=1

# #generate segments that connect different colors and record these?
colors = [0,1,2]
skipHoles = False
segments = set()
nodes    = set()
holes    = set()
#take care of boundary for domain first?
for j in range(0,ny):
    if aSlice[0,j] == colors[0]:
        nodes.add((0,j))  #bottom corner
        nodes.add((0,j+1))#top corner
        segments.add(((0,j),(0,j+1)))
    if aSlice[nx-1,j] == colors[0]:
        nodes.add((nx,j))  #bottom corner
        nodes.add((nx,j+1))#top corner
        segments.add(((nx,j),(nx,j+1)))
    #
#
for i in range(0,nx):
    if aSlice[i,0] == colors[0]:
        nodes.add((i,0))  #left corner
        nodes.add((i+1,0))#right corner
        segments.add(((i,0),(i+1,0)))
    if aSlice[i,ny-1] == colors[0]:
        nodes.add((i,ny))  #left corner
        nodes.add((i+1,ny))#right corner
        segments.add(((i,ny),(i+1,ny)))
    #
#
# nodes.add((0,0))
# nodes.add((0,ny))
# nodes.add((nx,ny))
# nodes.add((nx,0))
# segments.add(((0,0),(0,ny)))
# segments.add(((0,ny),(nx,ny)))
# segments.add(((nx,ny),(nx,0)))
# segments.add(((nx,0),(0,0)))

#x faces
for i in range(nx-1):
    for j in range(ny):
        if (aSlice[i,j] == colors[0] and aSlice[i+1,j] in colors[1:] or
            (aSlice[i,j] in colors[1:] and aSlice[i+1,j] == colors[0])) :
            nodes.add((i+1,j))  #bottom corner
            nodes.add((i+1,j+1))#top corner
            segments.add(((i+1,j),(i+1,j+1)))
            #holes.add((i+1.5,j))

        #
    #
#y faces
for j in range(ny-1):
    for i in range(nx):
        if (aSlice[i,j] == colors[0] and aSlice[i,j+1] in colors[1:] or
            (aSlice[i,j] in colors[1:] and aSlice[i,j+1] == colors[0])):
            nodes.add((i,j+1))  #left corner
            nodes.add((i+1,j+1))#right corner
            segments.add(((i,j+1),(i+1,j+1)))
            #holes.add((i,j+1.5))

    #
#
#
nxh=2
nyh=2
for j in range(nyh,ny-nyh):
    for i in range(nxh,nx-nxh):
        if aSlice[i,j] == colors[1]:
            found0 = False
            for jj in range(-nyh,nyh):
                for ii in range(-nxh,nxh):
                    if aSlice[i+ii,j+jj] == colors[0]:
                        found0 = True
            if not found0:
                holes.add((i,j))
        #
    #
##
sout = open('trislice.poly','w')
sout.write("#nnodes 2 nattributes nboundary markers\n")
sout.write("%d 2 0 0\n" % len(nodes))
nvert=0;
nodeid= {}
for n in nodes:
    sout.write("%d %g %g\n" % (nvert+base,n[0]*dx,n[1]*dy))
    nodeid[n] = nvert
    nvert+=1
#
sout.write("#nseg nboundary markers\n")
sout.write("%d 0\n" % len(segments))
nseg = 0
for s in segments:
    sout.write("%d %d %d\n" % (nseg+base,nodeid[s[0]]+base,nodeid[s[1]]+base))
    nseg+= 1
#

sout.write("#nholes \n")
if skipHoles == True:
    sout.write("%d \n" % 0)
else:
    sout.write("%d \n" % len(holes))
    nhole=0
    for h in holes:
        sout.write("%d %g %g\n" % (nhole+base,h[0]*dx,h[1]*dy))
        nhole+=1
    #
#
sout.write("\n")
sout.close()
