import math
from pyadh import Domain
#This is the channel from the side
#              top
#             ____________
#            |            |
# upstream   | O          |  downstream
#            |            |
#             ____________
#             bottom
#^
#y x>

#
#an example of a function defining the cross section of the obstacle
#
def circular_cross_section(center,radius,theta):
    return (radius*math.cos(theta)+center[0],
            radius*math.sin(theta)+center[1])

#
# a function to build this domain
#
def squareWithHole2d(height=1.0,
                     length=1.0,
                     radius=0.25,
                     center=[0.5,0.5],
                     n_points_on_hole=10,
                     cross_section=circular_cross_section):
    #enumerate the boundaries so we can supply integer flags for vertices and segments
    boundaries=['left',
                'right',
                'bottom',
                'top',
                'hole']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #now work around the domain from (0.0,0.0) going counterclockwise building vertices and segments (and flags)
    vertices = [[0.0,0.0],
                [length,0.0],
                [length,height],
                [0.0,height]]
    vertexFlags = [boundaryTags['left'],
                   boundaryTags['right'],
                   boundaryTags['right'],
                   boundaryTags['left']]
    segments = [[0,1],
                [1,2],
                [2,3],
                [3,0]]
    segmentFlags = [boundaryTags['bottom'],
                    boundaryTags['right'],
                    boundaryTags['top'],
                    boundaryTags['left']]
    #now do the vertices and segments on the squareWithHole
    theta=0.0
    pb  = cross_section(center,radius,theta)
    vertices.append([pb[0],pb[1]])
    vertexFlags.append(boundaryTags['hole'])
    for gb in range(1,n_points_on_hole):
        theta = float(gb)/float(n_points_on_hole)*2.0*math.pi
        pb  = cross_section(center,radius,theta)
        vertices.append([pb[0],pb[1]])
        vertexFlags.append(boundaryTags['hole'])
        nv = len(vertices)
        segments.append((nv-2,nv-1))
        segmentFlags.append(boundaryTags['hole'])
    segments.append((nv-1,nv-n_points_on_hole))
    segmentFlags.append(boundaryTags['hole'])
    #done with vertices and segments
    #holes
    holes=[center]
    #regions (just pick a point inside to indicate the flow region)
    regions=[[vertices[0][0]+length*1.0e-8,
              vertices[0][1]+length*1.0e-8]]
    regionFlags=[1]
    #construct domain object
    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  holes=holes,
                                                  regions=regions,
                                                  regionFlags=regionFlags)
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain =  squareWithHole2d(height=2.5,length=10.1,radius=0.1,cross_section=circular_cross_section,center=[0.9,1.3])
    domain.writeAsymptote("squareWithHole2d")
    domain.writePoly("squareWithHole2d")
    os.system("asy -V squareWithHole2d")
