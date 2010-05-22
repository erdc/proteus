#!/usr/bin/env python

import math
def genPoly(fileprefix,
            nx,
            ny,
            points_on_grain=50,
            points_on_boundary=50):
    """
    generate square lattice circle packing on unit square
    radius is determined by number of spheres
    returns boundary flags
    """
    #n_domain_vertices = 4
    #domain_vertices =[(0.0,0.0),(0.0,1.0),(1.0,1.0),(1.0,0.0)]
    n_domain_vertices = 4 + 4*points_on_boundary
    DX=1.0/float(points_on_boundary+1)
    DY=1.0/float(points_on_boundary+1)
    domain_vertices=[]
    for k in range(points_on_boundary+2):
        domain_vertices.append((0.0,k*DY))
    for k in range(1,points_on_boundary+2):
        domain_vertices.append((k*DX,1.0))
    for k in range(points_on_boundary,-1,-1):
        domain_vertices.append((1.0,k*DY))
    for k in range(points_on_boundary,0,-1):
        domain_vertices.append((k*DX,0.0))
    dx = 1.0/float(nx)
    dy = 1.0/float(ny)
    radius = 0.25*min(dx,dy)
    grain_centers  = []
    for i in range(ny):
        for j in range(nx):
            grain_centers.append((i*dy+0.5*dy,j*dx+0.5*dx))
    nvertices = len(domain_vertices) + len(grain_centers)*points_on_grain
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (nvertices,2,0,1))
    boundaries = ['left','right','front','back']
    #based on vertex order and way initial segment list goes round
    boundaryFlags={'left':1,'back':2,'right':3,'front':4} 
    #write vertices
    poly.write("#vertices \n")
    for v,p in enumerate(domain_vertices):#numbering is base 1 for triangle
        if (p[0] == 0.0):
            flag=boundaryFlags['left']
        elif (p[0] == 1.0):
            flag=boundaryFlags['right']
        elif (p[1] == 0.0):
            flag=boundaryFlags['front']
        elif (p[1] == 1.0):
            flag=boundaryFlags['back']
        else:
            exit
        poly.write('%d %18.12e %18.12e %d #exterior \n' % (v+1,p[0],p[1],flag))
    #write segments
    #left is X_minus, right is X_plus, front is Y_minus, back is Y_plus
    segments=[]
    segmentGroups=[]
    for sN in range(len(domain_vertices)-1):
        segments.append([sN,sN+1])
        if (domain_vertices[sN][0] == 0.0 and domain_vertices[sN+1][0] ==  0.0):
            segmentGroups.append(boundaryFlags['left'])
        elif (domain_vertices[sN][0] == 1.0 and domain_vertices[sN+1][0] ==  1.0):
            segmentGroups.append(boundaryFlags['right'])
        elif (domain_vertices[sN][1] == 0.0 and domain_vertices[sN+1][1] ==  0.0):
            segmentGroups.append(boundaryFlags['front'])
        elif (domain_vertices[sN][1] == 1.0 and domain_vertices[sN+1][1] ==  1.0):
            segmentGroups.append(boundaryFlags['back'])
        else:
            exit
    segments.append([len(domain_vertices)-1,0])
    if (domain_vertices[segments[-1][0]][0] == 0.0 and domain_vertices[segments[-1][1]][0] ==  0.0):
        segmentGroups.append(boundaryFlags['left'])
    if (domain_vertices[segments[-1][0]][0] == 1.0 and domain_vertices[segments[-1][1]][0] ==  1.0):
        segmentGroups.append(boundaryFlags['right'])
    if (domain_vertices[segments[-1][0]][1] == 0.0 and domain_vertices[segments[-1][1]][1] ==  0.0):
        segmentGroups.append(boundaryFlags['front'])
    if (domain_vertices[segments[-1][0]][1] == 1.0 and domain_vertices[segments[-1][1]][1] ==  1.0):
        segmentGroups.append(boundaryFlags['back'])
    #end exterior boundary segments
    vStart = len(domain_vertices)
    sStart = len(segments)
    for g,c in enumerate(grain_centers):
        for gb in range(points_on_grain):
            pb  = (radius*math.sin(float(gb)/float(points_on_grain)*2.0*math.pi),radius*math.cos(float(gb)/float(points_on_grain)*2.0*math.pi))
            poly.write('%d %18.12e %18.12e %d \n' % (vStart + gb+1,c[0]+pb[0],c[1]+pb[1],4+1+g))
        for gb in range(points_on_grain-1):
            segments.append([sStart+gb,sStart+gb+1])
            segmentGroups.append(4+g+1)
        segments.append([sStart+points_on_grain-1,sStart])
        segmentGroups.append(4+g+1)
        key = 'circle_%s' % g
        boundaryFlags[key] = 4+g+1
        vStart = vStart + points_on_grain
        sStart = sStart + points_on_grain
    poly.write('%d %d \n' % (len(segments),1))
    #print segments
    poly.write("#segments \n")
    for sN,s,sG in zip(range(len(segments)),segments,segmentGroups):
        poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,sG))
    poly.write('%d \n' % len(grain_centers))
    poly.write("#holes \n")
    for gN,g in enumerate(grain_centers):
        poly.write('%d %18.12e %18.12e \n' % (gN+1,g[0],g[1]))
    poly.write('#regions\n')
    poly.write('%d\n' % (0,))
    poly.close()

    return boundaryFlags

if __name__=='__main__':
    import os
    genPoly('pome_test',nx=2,ny=2,points_on_grain=10,points_on_boundary=10)
    os.system('triangle -YYpAq30Den pome_test.poly')
    os.system('showme pome_test.1.ele')
