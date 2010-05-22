import math
from pyadh import Domain

def pome3d(radius=1.0,
           nx=2,
           ny=2,
           nz=2,
           nSectors = 2,
           pad=1.0,
           pointFacets=False):
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

    #multiple of radius to offset from boundaries
    noff= pad*0.5+0.5
    #cube boundary vertices
    nN=0
    nF=0
    vertices=[]
    facets=[]
    vertices =[[0.0,0.0,0.0],                                           #0
               [0.0,(ny+pad)*diameter,0.0],                             #1 
               [(nx+pad)*diameter,(ny+pad)*diameter,0.0],               #2
               [(nx+pad)*diameter,0.0,0.0],                             #3
               [0.0,0.0,(nz+pad)*diameter],                             #4
               [0.0,(ny+pad)*diameter,(nz+pad)*diameter],               #5
               [(nx+pad)*diameter,(ny+pad)*diameter,(nz+pad)*diameter], #6
               [(nx+pad)*diameter,0.0,(nz+pad)*diameter]]               #7
    #domain size
    L = ((nx+pad)*diameter,
         (ny+pad)*diameter,
         (nz+pad)*diameter)

    nN = 8
    #cube boundary facets
    facets=[[[0,1,2,3]],#bottom 1
            [[0,1,5,4]],#left  2
            [[0,3,7,4]],#front   3
            [[1,2,6,5]],#back  4
            [[2,3,7,6]],#right   5
            [[4,5,6,7]]]#top    6
    facetFlags = [1,2,3,4,5,6]
    #
    regions=[[1.0e-8,1.0e-8,1.0e-8]]
    regionFlags=[1]
    #
    sphereTag = 10
    boundaries = ['bottom','left','front','back','right','top']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    boundaryTags['sphere']=sphereTag

    nF = 6
    #adding poles of grains to make sure those nodes are unique
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                grain_centers[(i,j,k)]=[k*diameter+noff*diameter,
                                        j*diameter+noff*diameter,
                                        i*diameter+noff*diameter]
                north_poles[(i,j,k)]=nN
                vertices.append([k*diameter+noff*diameter,
                                 j*diameter+noff*diameter,
                                 i*diameter+noff*diameter+radius])
                nN+=1
                back_poles[(i,j,k)]=nN
                vertices.append([k*diameter+noff*diameter,
                                 j*diameter+noff*diameter+radius,
                                 i*diameter+noff*diameter])
                nN+=1
                right_poles[(i,j,k)]=nN
                vertices.append([k*diameter+noff*diameter+radius,
                                 j*diameter+noff*diameter,
                                 i*diameter+noff*diameter])
                nN+=1
                if i > 0:
                    south_poles[(i,j,k)] = north_poles[(i-1,j,k)]
                else:
                    south_poles[(i,j,k)]=nN
                    vertices.append([k*diameter+noff*diameter,
                                     j*diameter+noff*diameter,
                                     i*diameter+noff*diameter-radius])
                    nN+=1
                if j > 0:
                    front_poles[(i,j,k)] = back_poles[(i,j-1,k)]
                else:
                    front_poles[(i,j,k)] = nN
                    vertices.append([k*diameter+noff*diameter,
                                     j*diameter+noff*diameter-radius,
                                     i*diameter+noff*diameter])
                    nN+=1
                if k > 0:
                    left_poles[(i,j,k)] = right_poles[(i,j,k-1)]
                else:
                    left_poles[(i,j,k)] = nN
                    vertices.append([k*diameter+noff*diameter-radius,
                                     j*diameter+noff*diameter,
                                     i*diameter+noff*diameter])
                    nN+=1
                if pointFacets:
                    if i==0:
                        facets[0].append([south_poles[(i,j,k)]])
                    if i==nz-1:
                        facets[5].append([north_poles[(i,j,k)]])
                    if j==0:
                        facets[2].append([front_poles[(i,j,k)]])
                    if j==ny-1:
                        facets[3].append([back_poles[(i,j,k)]])
                    if k==0:
                        facets[1].append([left_poles[(i,j,k)]])
                    if k==nx-1:
                        facets[4].append([right_poles[(i,j,k)]])
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
                            xc = radius*math.cos(thetac)*math.cos(eta)
                            yc = radius*math.sin(thetac)*math.cos(eta)
                            zc = radius*math.sin(eta)
                            #physical coordinates
                            vertices.append([xc+grain_centers[(i,j,k)][0],
                                             yc+grain_centers[(i,j,k)][1],
                                             zc+grain_centers[(i,j,k)][2]])
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
                            xc = radius*math.cos(thetac)*math.cos(eta)
                            yc = radius*math.sin(thetac)*math.cos(eta)
                            zc = -radius*math.sin(eta)
                            #print xc,yc,zc
                            #physical coordinates
                            vertices.append([xc+grain_centers[(i,j,k)][0],
                                             yc+grain_centers[(i,j,k)][1],
                                             zc+grain_centers[(i,j,k)][2]])
                            nN+=1
                for ii in range(2*nSectors):
                    for jj in range(2*nSectors):
                        if ((ii < nSectors and jj < nSectors) or
                            (ii>=nSectors and  jj>=nSectors)):
                            #top
                            facets.append([[top_nodes[(ii,jj)],top_nodes[(ii+1,jj)],top_nodes[(ii+1,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            facets.append([[top_nodes[(ii,jj)],top_nodes[(ii+1,jj+1)],top_nodes[(ii,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            nF+=2
                            #bottom
                            facets.append([[bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj)],bottom_nodes[(ii+1,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            facets.append([[bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj+1)],bottom_nodes[(ii,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            nF+=2
                        else:
                            #top
                            facets.append([[top_nodes[(ii,jj)],top_nodes[(ii+1,jj)],top_nodes[(ii,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            facets.append([[top_nodes[(ii+1,jj)],top_nodes[(ii+1,jj+1)],top_nodes[(ii,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            nF+=2
                            #bottom
                            facets.append([[bottom_nodes[(ii,jj)],bottom_nodes[(ii+1,jj)],bottom_nodes[(ii,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            facets.append([[bottom_nodes[(ii+1,jj)],bottom_nodes[(ii+1,jj+1)],bottom_nodes[(ii,jj+1)]]])
                            facetFlags.append(boundaryTags['sphere'])
                            nF+=2
    #print vertices
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 regions=regions,
                                                 regionFlags=regionFlags,
                                                 holes=grain_centers.values())
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain =  pome3d()
    domain.writeAsymptote("pome3d")
    domain.writePoly("pome3d")
    os.system("asy -V pome3d")
