#!/usr/bin/env python
import math

def genPoly(fileprefix,
            radius,
            nx,
            ny,
            nz,
            nxpad=0,
            nypad=0,
            nzpad=0,
            nSectors = 5):
    """
    regualar nx,ny,nz cube packing with no padding at boundary
    returns domain size and boundary flags
    nxpad,nypad,nzpad are the number of radii to add on either side
    of the spheres before putting the box boundary
    if pad = 0, then have to make sure we put in the poles on the boundary
    as part of the mesh
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

    
    #cube boundary vertices
    nN=0
    nF=0
    vertices=[]
    facets=[]
    #vertices= []
    vertices =[(0.0,0.0,0.0),                                                          #0
               (0.0,(ny+nypad)*diameter,0.0),                                       #1 
               ((nx+nxpad)*diameter,(ny+nypad)*diameter,0.0),                    #2
               ((nx+nxpad)*diameter,0.0,0.0),                                       #3
               (0.0,0.0,(nz+nzpad)*diameter),                                       #4
               (0.0,(ny+nypad)*diameter,(nz+nzpad)*diameter),                    #5
               ((nx+nxpad)*diameter,(ny+nypad)*diameter,(nz+nzpad)*diameter), #6
               ((nx+nxpad)*diameter,0.0,(nz+nzpad)*diameter)]                    #7
    nN = 8
    #reverse polish notation
    #indeces now go -1 ... nx since grain centers are 0,...nx-1 etc
    corner_vertices = {}
    corner_vertices[(-1,-1,-1)]        = 0
    corner_vertices[(-1,ny,-1)]      = 1
    corner_vertices[(-1,ny,nx)]    = 2
    corner_vertices[(-1,-1,nx)]      = 3
    corner_vertices[(nz,-1,-1)]      = 4
    corner_vertices[(nz,ny,-1)]    = 5
    corner_vertices[(nz,ny,nx)]    = 6
    corner_vertices[(nz,-1,nx)]    = 7

    #domain size
    L = ((nx + nxpad)*diameter, #get 2*nxpad*radius padding 
         (ny + nypad)*diameter,
         (nz + nzpad)*diameter)
    sphereTag = 10
    boundaries = ['bottom','top','left','right','front','back']#'front','left','right','back','top']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    boundaryTags['sphere']=sphereTag

    #cube boundary facets
    #facets=[(0,1,2,3),#bottom 1
    #        (0,1,5,4),#front  2
    #        (0,3,7,4),#left   3
    #        (1,2,6,5),#right  4
    #        (2,3,7,6),#back   5
    #        (4,5,6,7)]#top    6

    #adding poles of grains to make sure those nodes are unique
    #offset is base on pad
    nxoff = nxpad*radius; nyoff = nypad*radius; nzoff = nzpad*radius
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                grain_centers[(i,j,k)]=(k*diameter+nxoff+radius,
                                        j*diameter+nyoff+radius,
                                        i*diameter+nzoff+radius)
                north_poles[(i,j,k)]=nN
                vertices.append((k*diameter+nxoff+radius,
                                 j*diameter+nyoff+radius,
                                 i*diameter+nzoff+radius+radius))
                nN+=1
                back_poles[(i,j,k)]=nN
                vertices.append((k*diameter+nxoff+radius,
                                 j*diameter+nyoff+radius+radius,
                                 i*diameter+nzoff+radius))
                nN+=1
                right_poles[(i,j,k)]=nN
                vertices.append((k*diameter+nxoff+radius+radius,
                                 j*diameter+nyoff+radius,
                                 i*diameter+nzoff+radius))
                nN+=1
                if i > 0:
                    south_poles[(i,j,k)] = north_poles[(i-1,j,k)]
                else:
                    south_poles[(i,j,k)]=nN
                    vertices.append((k*diameter+nxoff+radius,
                                     j*diameter+nyoff+radius,
                                     i*diameter+nzoff+radius-radius))
                    nN+=1
                if j > 0:
                    front_poles[(i,j,k)] = back_poles[(i,j-1,k)]
                else:
                    front_poles[(i,j,k)] = nN
                    vertices.append((k*diameter+nxoff+radius,
                                     j*diameter+nyoff+radius-radius,
                                     i*diameter+nzoff+radius))
                    nN+=1
                if k > 0:
                    left_poles[(i,j,k)] = right_poles[(i,j,k-1)]
                else:
                    left_poles[(i,j,k)] = nN
                    vertices.append((k*diameter+nxoff+radius-radius,
                                     j*diameter+nyoff+radius,
                                     i*diameter+nzoff+radius))
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
    #

    boundaryFacesStart = nF
    #to start, add vertices even if there is padding
    #now add grid of vertices on boundary
    #initial ring on outer edges
    bp = set()
    for i in [-1,nz]:
        bp |= set([(i,j,k) for j in [-1,ny] for k in range(0,nx)])
        bp |= set([(i,j,k) for j in range(0,ny) for k in [-1,nx]])
    for j in [-1,ny]:
        bp |= set([(i,j,k) for i in [-1,nz] for k in range(0,nx)])
        bp |= set([(i,j,k) for i in range(0,nz) for k in [-1,nx]])
    for k in [-1,nx]:
        bp |= set([(i,j,k) for i in [-1,nz] for j in range(0,ny)])
        bp |= set([(i,j,k) for i in range(0,nz) for j in [-1,ny]])
    boundary_vertices = {}
    #mwf debug
    #import pdb
    #pdb.set_trace()
    for p in bp:
        if p[2] in [-1,nx]: #domain boundary
            x = (p[2]+1)*L[0]/(nx+1)
        else: #grain center 
            x = p[2]*diameter+radius+nxoff
        if p[1] in [ny,-1]: #domain boundary
            y = L[1]*(p[1]+1)/(ny+1)
        else:
            y = p[1]*diameter+radius+nyoff
        if p[0] in [nz,-1]: #domain boundary
            z = L[2]*(p[0]+1)/(nz+1)
        else:
            z =p[0]*diameter+radius+nzoff
        vertices.append((x,y,z))
        boundary_vertices[p] = nN
        nN += 1

    #rest of points correspond to grain centers projected onto face
    for i,z in zip([-1,nz],[0.,L[2]]):
        faceVertices = {}
        interior_points = [(i,j,k) for j in range(0,ny) for k in range(0,nx)]
        corner_points = [(i,j,k) for j in [-1,ny] for k in [-1,nx]]
        boundary_points = [(i,j,k) for j in [-1,ny] for k in range(0,nx)]
        boundary_points.extend([(i,j,k) for j in range(0,ny) for k in [-1,nx]])

        for p in corner_points:
            faceVertices[p] = corner_vertices[p]
        for p in boundary_points:
            faceVertices[p]=boundary_vertices[p]

        if nzpad == 0:
            for p in interior_points:
                if p[0] == -1:
                    faceVertices[p] = south_poles[(0,p[1],p[2])]
                else:
                    faceVertices[p] = north_poles[(nz-1,p[1],p[2])]
        else:
            for p in interior_points:
                x = p[2]*diameter + nxoff + radius
                y = p[1]*diameter + nyoff + radius
                vertices.append((x,y,z))
                faceVertices[p] = nN
                nN += 1
 
        for j in range(-1,ny):
            for k in range(-1,nx):
                facets.append((faceVertices[(i,j,k)],faceVertices[(i,j,k+1)],
                               faceVertices[(i,j+1,k+1)],faceVertices[(i,j+1,k)]))
                nF += 1
        #
    #
    for j,y in zip([-1,ny],[0.,L[1]]):
        faceVertices = {}
        interior_points = [(i,j,k) for i in range(0,nz) for k in range(0,nx)]
        corner_points = [(i,j,k) for i in [-1,nz] for k in [-1,nx]]
        boundary_points = [(i,j,k) for i in [-1,nz] for k in range(0,nx)]
        boundary_points.extend([(i,j,k) for i in range(0,nz) for k in [-1,nx]])
        for p in corner_points:
            faceVertices[p] = corner_vertices[p]
        for p in boundary_points:
            faceVertices[p]=boundary_vertices[p]
            #
        #
        #get interior nodes from poles
        if nypad == 0:
            for p in interior_points:
                if p[1] == -1:
                    faceVertices[p] = front_poles[(p[0],0,p[2])]
                else:
                    faceVertices[p] = back_poles[(p[0],ny-1,p[2])]
        else:
            for p in interior_points:
                x = p[2]*diameter + nxoff + radius
                z = p[0]*diameter + nzoff + radius
                vertices.append((x,y,z))
                faceVertices[p] = nN
                nN += 1
            
        for i in range(-1,nz):
            for k in range(-1,nx):
                facets.append((faceVertices[(i,j,k)],faceVertices[(i,j,k+1)],
                               faceVertices[(i+1,j,k+1)],faceVertices[(i+1,j,k)]))
                nF += 1
        #
    #
    for k,x in zip([-1,nx],[0.,L[0]]):
        faceVertices = {}
        interior_points = [(i,j,k) for i in range(0,nz) for j in range(0,ny)]
        corner_points = [(i,j,k) for i in [-1,nz] for j in [-1,ny]]
        boundary_points = [(i,j,k) for i in [-1,nz] for j in range(0,ny)]
        boundary_points.extend([(i,j,k) for i in range(0,nz) for j in [-1,ny]])
        for p in corner_points:
            faceVertices[p] = corner_vertices[p]
        for p in boundary_points:
            faceVertices[p]=boundary_vertices[p]
            #
        #
        #get interior nodes from poles
        if nxpad == 0:
            for p in interior_points:
                if p[2] == -1:
                    faceVertices[p] = left_poles[(p[0],p[1],0)]
                else:
                    faceVertices[p] = right_poles[(p[0],p[1],nx-1)]
        else:
           for p in interior_points:
                z = p[0]*diameter + nzoff + radius
                y = p[1]*diameter + nyoff + radius
                vertices.append((x,y,z))
                faceVertices[p] = nN
                nN += 1
            
        for i in range(-1,nz):
            for j in range(-1,ny):
                facets.append((faceVertices[(i,j,k)],faceVertices[(i,j+1,k)],
                               faceVertices[(i+1,j+1,k)],faceVertices[(i+1,j,k)]))
                nF += 1
        #
    #


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
    for fN,f in enumerate(facets[:boundaryFacesStart]):
        poly.write('%d %d %d \n' % (1,0,sphereTag)) #1 polygon, 0 holes, 1 on boudary
        poly.write('%d %d %d %d \n' % (3,f[0]+1,f[1]+1,f[2]+1)) #facet number and corners  of facet
        #no holes in facet

    offset=boundaryFacesStart
    for i in range(2):
        nfaces = (nx+1)*(ny+1)
        for fN,f in enumerate(facets[offset:offset+nfaces]):
            poly.write('%d %d %d \n' % (1,0,i+1)) #1 polygon, 0 holes, 1 on boudary
            poly.write('%d %d %d %d %d\n' % (4,f[0]+1,f[1]+1,f[2]+1,f[3]+1)) #facet number and corners  of facet
        #no holes in facet
        offset += nfaces
    for j in range(2):
        nfaces = (nx+1)*(nz+1)
        for fN,f in enumerate(facets[offset:offset+nfaces]):
            poly.write('%d %d %d \n' % (1,0,3+j)) #1 polygon, 0 holes, 1 on boudary
            poly.write('%d %d %d %d %d\n' % (4,f[0]+1,f[1]+1,f[2]+1,f[3]+1)) #facet number and corners  of facet
        #no holes in facet
        offset += nfaces
    for k in range(2):
        nfaces = (ny+1)*(nz+1)
        for fN,f in enumerate(facets[offset:offset+nfaces]):
            poly.write('%d %d %d \n' % (1,0,5+k)) #1 polygon, 0 holes, 1 on boudary
            poly.write('%d %d %d %d %d\n' % (4,f[0]+1,f[1]+1,f[2]+1,f[3]+1)) #facet number and corners  of facet
        #no holes in facet
        offset += nfaces
    assert offset == len(facets)

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
    radius = 1.; nx = 4; ny = 4; nz = 4; nxpad = 3; nypad = 3; nzpad = 3;
    L,bt = genPoly("pome_box",
                   radius,
                   nx,
                   ny,
                   nz,
                   nxpad,
                   nypad,
                   nzpad)

    print "radius= %s nx=%s ny=%s nz=%s nxpad=%s nypad=%s nzpad=%s L=%s " % (radius,
                                                                             nx,ny,nz,
                                                                             nxpad,nypad,nzpad,L)
