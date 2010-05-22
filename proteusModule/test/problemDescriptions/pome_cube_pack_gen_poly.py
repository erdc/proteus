#!/usr/bin/env python
import math

def genPoly(fileprefix,
            radius,
            nx,
            ny,
            nz,
            nSectors = 5):
    """
    regualar nx,ny,nz cube packing with padding to boundary
    returns domain size and boundary flags
    """
    n_domain_vertices = 8
    diameter=2.0*radius
    grain_centers  = {}
    north_poles = {}
    right_poles = {}
    back_poles = {}
    south_poles = {}
    left_poles = {}
    front_poles = {}

    pad= 3
    noff= pad*0.5+0.5
    #cube boundary vertices
    nN=0
    nF=0
    vertices=[]
    facets=[]
    vertices =[(0.0,0.0,0.0),                                     #0
               (0.0,(ny+pad)*diameter,0.0),                         #1 
               ((nx+pad)*diameter,(ny+pad)*diameter,0.0),             #2
               ((nx+pad)*diameter,0.0,0.0),                         #3
               (0.0,0.0,(nz+pad)*diameter),                         #4
               (0.0,(ny+pad)*diameter,(nz+pad)*diameter),             #5
               ((nx+pad)*diameter,(ny+pad)*diameter,(nz+pad)*diameter), #6
               ((nx+pad)*diameter,0.0,(nz+pad)*diameter)]             #7
    #domain size
    L = ((nx+pad)*diameter,
         (ny+pad)*diameter,
         (nz+pad)*diameter)

    nN = 8
    #cube boundary facets
    facets=[(0,1,2,3),#bottom 1
            (0,1,5,4),#front  2
            (0,3,7,4),#left   3
            (1,2,6,5),#right  4
            (2,3,7,6),#back   5
            (4,5,6,7)]#top    6
    sphereTag = 10
    boundaries = ['bottom','front','left','right','back','top']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    boundaryTags['sphere']=sphereTag

    nF = 6
    #adding poles of grains to make sure those nodes are unique
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                grain_centers[(i,j,k)]=(k*diameter+noff*diameter,
                                        j*diameter+noff*diameter,
                                        i*diameter+noff*diameter)
                north_poles[(i,j,k)]=nN
                vertices.append((k*diameter+noff*diameter,
                                 j*diameter+noff*diameter,
                                 i*diameter+noff*diameter+radius))
                nN+=1
                back_poles[(i,j,k)]=nN
                vertices.append((k*diameter+noff*diameter,
                                 j*diameter+noff*diameter+radius,
                                 i*diameter+noff*diameter))
                nN+=1
                right_poles[(i,j,k)]=nN
                vertices.append((k*diameter+noff*diameter+radius,
                                 j*diameter+noff*diameter,
                                 i*diameter+noff*diameter))
                nN+=1
                if i > 0:
                    south_poles[(i,j,k)] = north_poles[(i-1,j,k)]
                else:
                    south_poles[(i,j,k)]=nN
                    vertices.append((k*diameter+noff*diameter,
                                     j*diameter+noff*diameter,
                                     i*diameter+noff*diameter-radius))
                    nN+=1
                if j > 0:
                    front_poles[(i,j,k)] = back_poles[(i,j-1,k)]
                else:
                    front_poles[(i,j,k)] = nN
                    vertices.append((k*diameter+noff*diameter,
                                     j*diameter+noff*diameter-radius,
                                     i*diameter+noff*diameter))
                    nN+=1
                if k > 0:
                    left_poles[(i,j,k)] = right_poles[(i,j,k-1)]
                else:
                    left_poles[(i,j,k)] = nN
                    vertices.append((k*diameter+noff*diameter-radius,
                                     j*diameter+noff*diameter,
                                     i*diameter+noff*diameter))
                    nN+=1

    hxi = radius/(math.sqrt(2.0)*float(nSectors));
    heta = radius/(math.sqrt(2.0)*float(nSectors));
    #now loop over grains
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                #top half  sphere nodes
                top_nodes = {}
                for ii in range(2*nSectors+1):
                    for jj in range(2*nSectors+1):
                        if (ii,jj) == (0,0):
                            top_nodes[(ii,jj)] = left_poles[(i,j,k)]
                        elif (ii,jj) == (2*nSectors,0):
                            top_nodes[(ii,jj)] = back_poles[(i,j,k)]
                        elif (ii,jj) == (2*nSectors,2*nSectors):
                            top_nodes[(ii,jj)] = right_poles[(i,j,k)]
                        elif (ii,jj) == (0,2*nSectors):
                            top_nodes[(ii,jj)] = front_poles[(i,j,k)]
                        elif (ii,jj) == (nSectors,nSectors):
                            top_nodes[(ii,jj)] = north_poles[(i,j,k)]
                        else:
                            top_nodes[(ii,jj)] = nN
                            #rotate  corners of ref square to line up with poles
                            x0s = (jj-nSectors)*hxi
                            y0s = (ii-nSectors)*heta
                            r0s = math.sqrt(x0s**2 + y0s**2)
                            theta0s = math.atan2(y0s,x0s)
                            theta1s = theta0s - math.pi/4.0
                            r1s = r0s
                            x1s = r1s*math.cos(theta1s)
                            y1s = r1s*math.sin(theta1s)
                            #do each quadrant
                            if x1s >= 0.0  and y1s >=0.0:
                                rc = x1s + y1s
                                thetac = theta1s
                            elif x1s < 0.0 and y1s >=0.0:
                                rc = y1s - x1s
                                thetac = theta1s
                            elif x1s <= 0.0 and y1s < 0.0:
                                rc = -(x1s+y1s)
                                thetac = theta1s
                            else:
                                rc = x1s - y1s
                                thetac = theta1s
                            eta = 0.5*math.pi*(radius-rc)/radius
    #                         xc = rc*math.cos(thetac)
    #                         yc = rc*math.sin(thetac)
    #                         zc = math.sqrt(math.fabs(radius**2 - xc**2 - yc**2))
                            xc = radius*math.cos(thetac)*math.cos(eta)
                            yc = radius*math.sin(thetac)*math.cos(eta)
                            zc = radius*math.sin(eta)
                            #zc = math.sqrt(math.fabs(radius**2 - xc**2 - yc**2))
                            #print xc,yc,zc,rc,radius**2 - xc**2 - yc**2
                            #physical coordinates
                            vertices.append((xc+grain_centers[(i,j,k)][0],
                                             yc+grain_centers[(i,j,k)][1],
                                             zc+grain_centers[(i,j,k)][2]))
                            nN+=1
                #bottom half sphere nodes
                bottom_nodes = {}
                for ii in range(2*nSectors+1):
                    for jj in range(2*nSectors+1):
                        if (ii,jj) == (0,0):
                            bottom_nodes[(ii,jj)] = left_poles[(i,j,k)]
                        elif (ii,jj) == (2*nSectors,0):
                            bottom_nodes[(ii,jj)] = back_poles[(i,j,k)]
                        elif (ii,jj) == (2*nSectors,2*nSectors):
                            bottom_nodes[(ii,jj)] = right_poles[(i,j,k)]
                        elif (ii,jj) == (0,2*nSectors):
                            bottom_nodes[(ii,jj)] = front_poles[(i,j,k)]
                        elif (ii,jj) == (nSectors,nSectors):
                            bottom_nodes[(ii,jj)] = south_poles[(i,j,k)]
                        elif (ii in [0,2*nSectors] or
                              jj in [0,2*nSectors]):#equator
                            bottom_nodes[(ii,jj)] = top_nodes[(ii,jj)]
                        else:
                            bottom_nodes[(ii,jj)] = nN
                            #rotate  corners of ref square to line up with poles
                            x0s = (jj-nSectors)*hxi
                            y0s = (ii-nSectors)*heta
                            r0s = math.sqrt(x0s**2 + y0s**2)
                            theta0s = math.atan2(y0s,x0s)
                            theta1s = theta0s - math.pi/4.0
                            r1s = r0s
                            x1s = r1s*math.cos(theta1s)
                            y1s = r1s*math.sin(theta1s)
                            #do each quadrant
                            if x1s >= 0.0  and y1s >=0.0:
                                rc = x1s + y1s
                                thetac = theta1s
                            elif x1s < 0.0 and y1s >=0.0:
                                rc = y1s - x1s
                                thetac = theta1s
                            elif x1s <= 0.0 and y1s < 0.0:
                                rc = -(x1s+y1s)
                                thetac = theta1s
                            else:
                                rc = x1s - y1s
                                thetac = theta1s
                            eta = 0.5*math.pi*(radius-rc)/radius
    #                         xc = rc*math.cos(thetac)
    #                         yc = rc*math.sin(thetac)
    #                         zc = -math.sqrt(math.fabs(radius**2 - xc**2 - yc**2))
                            xc = radius*math.cos(thetac)*math.cos(eta)
                            yc = radius*math.sin(thetac)*math.cos(eta)
                            zc = -radius*math.sin(eta)
                            #print xc,yc,zc
                            #physical coordinates
                            vertices.append((xc+grain_centers[(i,j,k)][0],
                                             yc+grain_centers[(i,j,k)][1],
                                             zc+grain_centers[(i,j,k)][2]))
                            nN+=1
    #             for ii in range(2*nSectors+1):
    #                 line=""
    #                 for jj in range(2*nSectors+1):
    #                     line+=`top_nodes[(ii,jj)]`+'\t'
    #                 print line
    #             print "-------------------"
    #             for ii in range(2*nSectors+1):
    #                 line=""
    #                 for jj in range(2*nSectors+1):
    #                     line+=`bottom_nodes[(ii,jj)]`+'\t'
    #                 print line
                for ii in range(2*nSectors):
                    for jj in range(2*nSectors):
                        if ((ii < nSectors and jj < nSectors) or
                            (ii>=nSectors and  jj>=nSectors)):
                            #top
                            facets.append([top_nodes[(ii,jj)],top_nodes[(ii+1,jj)],top_nodes[(ii+1,jj+1)]])
                            facets.append([top_nodes[(ii,jj)],top_nodes[(ii+1,jj+1)],top_nodes[(ii,jj+1)]])
                            nF+=2
                            #bottom
                            facets.append([bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj)],bottom_nodes[(ii+1,jj+1)]])
                            facets.append([bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj+1)],bottom_nodes[(ii,jj+1)]])
                            nF+=2
                        else:
                            #top
                            facets.append([top_nodes[(ii,jj)],top_nodes[(ii+1,jj)],top_nodes[(ii,jj+1)]])
                            facets.append([top_nodes[(ii+1,jj)],top_nodes[(ii+1,jj+1)],top_nodes[(ii,jj+1)]])
                            nF+=2
                            #bottom
                            facets.append([bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj)],bottom_nodes[(ii,jj+1)]])
                            facets.append([bottom_nodes[(ii+1,jj)],bottom_nodes[(ii+1,jj+1)],bottom_nodes[(ii,jj+1)]])
                            nF+=2

    poly = open(fileprefix+'.poly','w')
    poly.write("#vertices \n")
    assert(nN == len(vertices)),"nN != nVertices"
    poly.write('%d %d %d %d \n' % (nN,3,0,1))
    for v,p in enumerate(vertices):
        poly.write('%d %18.12e %18.12e %18.12e %d \n' % (v+1,p[0],p[1],p[2],1))
    #write facets
    poly.write("#facets \n")
    assert(nF == len(facets)),"nF != len(facets)"
    poly.write('%d %d \n' % (nF,1))
    for fN,f in enumerate(facets[:6]):
        poly.write('%d %d %d \n' % (1,0,fN+1)) #1 polygon, 0 holes, 1 on boundary, check consistent with boundaryTags
        poly.write('%d %d %d %d %d\n' % (4,f[0]+1,f[1]+1,f[2]+1,f[3]+1)) #facet number and corners  of facet
        #no holes in facet
    for fN,f in enumerate(facets[6:]):
        poly.write('%d %d %d \n' % (1,0,sphereTag)) #1 polygon, 0 holes, 10 on boudary
        poly.write('%d %d %d %d \n' % (3,f[0]+1,f[1]+1,f[2]+1)) #facet number and corners  of facet
        #no holes in facet
    #write holes
    poly.write("#holes \n")
    poly.write('%d \n' % len(grain_centers))
    # print grain_centers
    for gN,g in enumerate(grain_centers.values()):
        poly.write('%d %18.12e %18.12e %18.12e\n' % (gN+1,g[0],g[1],g[2]))
    #poly.write("#holes \n")
    #poly.write('%d \n' % (0,))
    poly.close()

    return L,boundaryTags
if __name__ == '__main__':
    radius = 1.; nx = 4; ny = 4; nz = 4;
    L,bt = genPoly("pome_gen",
                   radius,
                   nx,
                   ny,
                   nz)

    print "radius= %s nx=%s ny=%s nz=%s L=%s " % (radius,nx,ny,nz,L)

