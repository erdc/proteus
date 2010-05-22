import math
from pyadh import Domain

def wigley3d(fileprefix,
             height=10.0,
             length=100.0,
             width=25.0,
             hull_beam = 10.0,
             hull_draft = 8.0,
             hull_length=  50.0,
             hull_center=(100.0/2.0,25.0/2.0,10.0/2.0),
             n_points_draft=3,
             n_points_length=3):
    boundaries=['left',
                'right',
                'bottom',
                'top',
                'front',
                'back',
                'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [[0.0,0.0,0.0],#0
                [length,0.0,0.0],#1
                [length,0.0,height],#2
                [0.0,0.0,height],#3
                [0.0,width,0.0],#4
                [length,width,0.0],#5
                [length,width,height],#6
                [0.0,width,height]]#7
    vertexFlags = [boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top'],
                   boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top']]
    facets =[[[0,1,2,3]],#front
             [[4,5,6,7]],#back
             [[0,1,5,4]],#bottom
             [[2,3,7,6]],#top
             [[0,3,7,4]],#left
             [[1,2,6,5]]]#right
    facetFlags=[boundaryTags['front'],
                boundaryTags['back'],
                boundaryTags['bottom'],
                boundaryTags['top'],
                boundaryTags['left'],
                boundaryTags['right']]
    holes=[]
    #wigley hull parameters
    dx = hull_length/float(n_points_length-1)
    dz = hull_draft/float(n_points_draft-1)
    #grid on right half of hull
    for i in range(n_points_length):
        for j in range(n_points_draft):
            x = i*dx 
            z = j*dz
            y = 0.5*hull_beam*(1.0 - (2.0*(x-0.5*hull_length)/hull_length)**2) * (1.0 - ((hull_draft-z)/hull_draft)**2)
            vertices.append([x+hull_center[0]-0.5*hull_length,
                             y+hull_center[1],
                             z+hull_center[2]-0.5*hull_draft])
            vertexFlags.append(boundaryTags['obstacle'])
    def vN_right(i,j):
        return 8 + i*n_points_draft+j
    for i in range(n_points_length-1):
        for j in range(n_points_draft-1):
            if i < n_points_length/2:
                facets.append([[vN_right(i,j),vN_right(i+1,j+1),vN_right(i+1,j)]])
                facetFlags.append(boundaryTags['obstacle'])
                facets.append([[vN_right(i,j),vN_right(i,j+1),vN_right(i+1,j+1)]])
                facetFlags.append(boundaryTags['obstacle'])
            else:
                facets.append([[vN_right(i,j),vN_right(i,j+1),vN_right(i+1,j)]])
                facetFlags.append(boundaryTags['obstacle'])
                facets.append([[vN_right(i,j+1),vN_right(i+1,j+1),vN_right(i+1,j)]])
                facetFlags.append(boundaryTags['obstacle'])                
    #grid on left half of hull
    for i in range(1,n_points_length-1):
        for j in range(1,n_points_draft):
            x = i*dx
            z = j*dz
            y = 0.5*hull_beam*(1.0 - (2.0*(x-0.5*hull_length)/hull_length)**2)*(1.0 - ((hull_draft-z)/hull_draft)**2)
            vertices.append([x+hull_center[0]-0.5*hull_length,
                             hull_center[1]-y,
                             z+hull_center[2]-0.5*hull_draft])
            vertexFlags.append(boundaryTags['obstacle'])
    def vN_left(i,j):
        if i== 0 or j==0:
            return vN_right(i,j)
        if i == (n_points_length-1):# or j==(n_points_draft-1):
            return vN_right(i,j)
        else:
            return 8 + n_points_length*n_points_draft+(i-1)*(n_points_draft-1)+j-1
    for i in range(n_points_length-1):
        for j in range(n_points_draft-1):
            if i < n_points_length/2:
                facets.append([[vN_left(i,j),vN_left(i+1,j+1),vN_left(i+1,j)]])
                facetFlags.append(boundaryTags['obstacle'])
                facets.append([[vN_left(i,j),vN_left(i,j+1),vN_left(i+1,j+1)]])
                facetFlags.append(boundaryTags['obstacle'])
            else:
                facets.append([[vN_left(i,j),vN_left(i,j+1),vN_left(i+1,j)]])
                facetFlags.append(boundaryTags['obstacle'])
                facets.append([[vN_left(i,j+1),vN_left(i+1,j+1),vN_left(i+1,j)]])
                facetFlags.append(boundaryTags['obstacle'])                
    topFacet=[]
    for i in range(n_points_length):
        topFacet.append(vN_right(i,n_points_draft-1))
    for i in range(n_points_length-1,0,-1):
        topFacet.append(vN_left(i,n_points_draft-1))
    facets.append([topFacet])
    facetFlags.append(boundaryTags['obstacle'])
    #for v in vertices: print v
    #for f in facets: print f
    holes.append(hull_center)
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 vertexFlags=vertexFlags,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 holes=holes)
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain = wigley3d(fileprefix="wigley3dDomainTest")
    domain.writePoly("wigley3dDomainTest")
    domain.writePLY("wigley3dDomainTest")
    domain.writeAsymptote("wigley3dDomainTest")
    os.system("asy -V wigley3dDomainTest.asy")
