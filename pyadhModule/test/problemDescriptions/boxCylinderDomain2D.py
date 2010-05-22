import math
from pyadh import Domain
#This is the channel from the side
#              top
#             _____________
#            |             |
# upstream   | O           |  downstream
#            |     _       |
#             ____| |_______
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
def boxCylinder2D(height=1.0,
                  length=1.0,
                  radius=0.25,
                  center=[0.5,0.5],
                  boxStart=0.5,
                  boxWidth=0.1,
                  boxHeight=0.1,
                  n_points_on_obstacle=10,
                  cross_section=circular_cross_section):
    #enumerate the boundaries so we can supply integer flags for vertices and segments
    boundaries=['left',
                'right',
                'bottom',
                'top',
                'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #now work around the domain from (0.0,0.0) going counterclockwise building vertices and segments (and flags)
    vertices = [[0.0,0.0],
                [length,0.0],
                [length,height],
                [0.0,height],
                [boxStart,0.0],
                [boxStart,boxHeight],
                [boxStart+boxWidth,boxHeight],
                [boxStart+boxWidth,0.0]]
    vertexFlags = [boundaryTags['bottom'], #0
                   boundaryTags['bottom'], #1
                   boundaryTags['top'],#2
                   boundaryTags['top'],#3
                   boundaryTags['bottom'],#4
                   boundaryTags['bottom'],#5
                   boundaryTags['bottom'],#6
                   boundaryTags['bottom']]#7
    segments = [[0,4],
                [4,5],
                [5,6],
                [6,7],
                [7,1],
                [1,2],
                [2,3],
                [3,0]]
    segmentFlags = [boundaryTags['bottom'],
                    boundaryTags['bottom'],
                    boundaryTags['bottom'],
                    boundaryTags['bottom'],
                    boundaryTags['bottom'],
                    boundaryTags['right'],
                    boundaryTags['top'],
                    boundaryTags['left']]
    #now do the vertices and segments on the cylinder
    theta=0.0
    pb  = cross_section(center,radius,theta)
    vertices.append([pb[0],pb[1]])
    vertexFlags.append(boundaryTags['obstacle'])
    for gb in range(1,n_points_on_obstacle):
        theta = float(gb)/float(n_points_on_obstacle)*2.0*math.pi
        pb  = cross_section(center,radius,theta)
        vertices.append([pb[0],pb[1]])
        vertexFlags.append(boundaryTags['obstacle'])
        nv = len(vertices)
        segments.append((nv-2,nv-1))
        segmentFlags.append(boundaryTags['obstacle'])
    segments.append((nv-1,nv-n_points_on_obstacle))
    segmentFlags.append(boundaryTags['obstacle'])
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
    domain =  boxCylinder2D(height=2.0,length=10.0,radius=0.5,cross_section=circular_cross_section,center=[4.0,1.0],boxStart=2.0,boxWidth=0.5,boxHeight=0.5)
    domain.writeAsymptote("boxCylinder2D")
    domain.writePoly("boxCylinder2D")
    os.system("asy -V boxCylinder2D")
