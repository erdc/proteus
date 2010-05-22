import math
from pyadh import Domain
#This is the channel from the side, rotated clockwise 90 degrees 
#              up
#             ____
#            |    |
#bottom      | O  |  top
#            |    |
#             ____
#             down
#y ->
#x
#|
#V

def circular_cross_section(center,radius,theta):
    return (radius*math.sin(theta)+center[0],
            radius*math.cos(theta)+center[2])

def cylinder3d(fileprefix,
               height=1.0,
               length=1.0,
               width=1.0,
               radius=0.1,
               center=(0.5,0.5,0.5),
               n_points_on_obstacle=2*21-2,
               cross_section=circular_cross_section,
               thetaOffset=0.0):
    boundaries=['upstream',
                'downstream',
                'bottom',
                'top',
                'front',
                'back',
                'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertexKeys = ['upstream_bottom',
                  'downstream_bottom',
                  'downstream_top',
                  'upstream_top']
    vertices = [[0.0,0.0],
                [length,0.0],
                [length,height],
                [0.0,height]]
    vertexFlags = [boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top']]
    nv = len(vertices)
    #cylinder
    theta=thetaOffset
    pb  = cross_section(center,radius,theta)
    vertices.append([pb[0],pb[1]])
    vertexKeys.append('obstacle_'+`0`)
    vertexFlags.append(boundaryTags['obstacle'])
    for gb in range(1,n_points_on_obstacle):
        theta = float(gb)/float(n_points_on_obstacle)*2.0*math.pi+thetaOffset
        pb  = cross_section(center,radius,theta)
        vertexKeys.append('obstacle_'+`gb`)
        vertices.append([pb[0],pb[1]])
        vertexFlags.append(boundaryTags['obstacle'])
    #now need to convert rep to 3D
    vertices3dDict={}
    vertices3d=[]
    vertexFlags3d=[]
    facets3d=[]
    facetFlags3d=[]
    facetHoles3d=[]
    front_cylinder=[]
    back_cylinder=[]
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_front'] = vN
        vertices3d.append([v[0],0.0,v[1]])
        vertexFlags3d.append(vertexFlags[vN])
        if 'obstacle' in vertexKeys[vN]:
            front_cylinder.append(vN)
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_back']=vN+len(vertices)
        vertices3d.append([v[0],width,v[1]])
        vertexFlags3d.append(vertexFlags[vN])
        if 'obstacle' in vertexKeys[vN]:
            back_cylinder.append(vN+len(vertices))
    #upstream
    facets3d.append([[vertices3dDict['upstream_bottom_front'],
                      vertices3dDict['upstream_bottom_back'],
                      vertices3dDict['upstream_top_back'],
                      vertices3dDict['upstream_top_front']]])
    facetFlags3d.append(boundaryTags['upstream'])
    facetHoles3d.append([])
    #downstream
    facets3d.append([[vertices3dDict['downstream_bottom_front'],
                     vertices3dDict['downstream_bottom_back'],
                     vertices3dDict['downstream_top_back'],
                     vertices3dDict['downstream_top_front']]])
    facetFlags3d.append(boundaryTags['downstream'])
    facetHoles3d.append([])
    #top
    facets3d.append([[vertices3dDict['upstream_top_front'],
                     vertices3dDict['downstream_top_front'],
                     vertices3dDict['downstream_top_back'],
                     vertices3dDict['upstream_top_back']]])
    facetFlags3d.append(boundaryTags['top'])
    facetHoles3d.append([])
    #bottom
    facets3d.append([[vertices3dDict['upstream_bottom_front'],
                     vertices3dDict['downstream_bottom_front'],
                     vertices3dDict['downstream_bottom_back'],
                     vertices3dDict['upstream_bottom_back']]])
    facetFlags3d.append(boundaryTags['bottom'])
    facetHoles3d.append([])
    #front
    facets3d.append([[vertices3dDict['upstream_bottom_front'],
                      vertices3dDict['downstream_bottom_front'],
                      vertices3dDict['downstream_top_front'],
                      vertices3dDict['upstream_top_front']],
                     front_cylinder])
    facetFlags3d.append(boundaryTags['front'])
    facetHoles3d.append([(center[0],0.0,center[2])])
    #back
    facets3d.append([[vertices3dDict['upstream_bottom_back'],
                     vertices3dDict['downstream_bottom_back'],
                     vertices3dDict['downstream_top_back'],
                     vertices3dDict['upstream_top_back']],
                     back_cylinder])
    facetFlags3d.append(boundaryTags['back'])
    facetHoles3d.append([(center[0],width,center[2])])
    #cylinder
    for fN in range(n_points_on_obstacle-1):
        facets3d.append([[front_cylinder[fN],
                         back_cylinder[fN],
                         back_cylinder[fN+1],
                         front_cylinder[fN+1]]])
        facetFlags3d.append(boundaryTags['obstacle'])
        facetHoles3d.append([])
    facets3d.append([[front_cylinder[-1],
                      back_cylinder[-1],
                      back_cylinder[0],
                      front_cylinder[0]]])
    facetFlags3d.append(boundaryTags['obstacle'])
    facetHoles3d.append([])
    print vertices3d,vertexFlags3d,facets3d,facetFlags3d,facetHoles3d
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices3d,
                                                 vertexFlags=vertexFlags3d,
                                                 facets=facets3d,
                                                 facetFlags=facetFlags3d,
                                                 facetHoles=facetHoles3d)
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain = cylinder3d(fileprefix="cylinder3dDomainTest",
                        height=1.0,
                        length=1.0,
                        width=1.0,
                        radius=0.1,
                        center=(0.5,0.5,0.5),
                        n_points_on_obstacle=2*5-2,
                        cross_section=circular_cross_section,
                        thetaOffset=0.0)
    domain.writeAsymptote("cylinder3dDomainTest")
    os.system("asy -V cylinder3dDomainTest")

