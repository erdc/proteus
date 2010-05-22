import math
from pyadh import Domain

def circular_cross_section(center,radius,theta):
    return (radius*math.sin(theta)+center[0],
            radius*math.cos(theta)+center[2])

def collector3d(fileprefix,
                length=1.0,
                radius=0.1,
                center=(0.5,0.5,0.5),
                n_points_on_screen=2*21-2,
                cross_section=circular_cross_section,
                thetaOffset=0.0):
    boundaries=['upstream',
                'downstream',
                'screen']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertexKeys = []
    vertices = []
    vertexFlags = []
    nv =0
    #cylinder
    theta=thetaOffset
    pb  = cross_section(center,radius,theta)
    vertices.append([pb[0],pb[1]])
    vertexKeys.append('screen_'+`0`)
    vertexFlags.append(boundaryTags['screen'])
    for gb in range(1,n_points_on_screen):
        theta = float(gb)/float(n_points_on_screen)*2.0*math.pi+thetaOffset
        pb  = cross_section(center,radius,theta)
        vertexKeys.append('screen_'+`gb`)
        vertices.append([pb[0],pb[1]])
        vertexFlags.append(boundaryTags['screen'])
    #now need to convert rep to 3D
    vertices3dDict={}
    vertices3d=[]
    vertexFlags3d=[]
    facets3d=[]
    facetFlags3d=[]
    facetHoles3d=[]
    upstream_cylinder=[]
    downstream_cylinder=[]
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_upstream'] = vN
        vertices3d.append([v[0],0.0,v[1]])
        vertexFlags3d.append(vertexFlags[vN])
        if 'screen' in vertexKeys[vN]:
            upstream_cylinder.append(vN)
    for vN,v in enumerate(vertices):
        vertices3dDict[vertexKeys[vN]+'_downstream']=vN+len(vertices)
        vertices3d.append([v[0],length,v[1]])
        vertexFlags3d.append(vertexFlags[vN])
        if 'screen' in vertexKeys[vN]:
            downstream_cylinder.append(vN+len(vertices))
    #upstream
    facets3d.append([upstream_cylinder])
    facetFlags3d.append(boundaryTags['upstream'])
    #downstream
    facets3d.append([downstream_cylinder])
    facetFlags3d.append(boundaryTags['downstream'])
    #cylinder
    for fN in range(n_points_on_screen-1):
        facets3d.append([[upstream_cylinder[fN],
                         downstream_cylinder[fN],
                         downstream_cylinder[fN+1],
                         upstream_cylinder[fN+1]]])
        facetFlags3d.append(boundaryTags['screen'])
        facetHoles3d.append([])
    facets3d.append([[upstream_cylinder[-1],
                      downstream_cylinder[-1],
                      downstream_cylinder[0],
                      upstream_cylinder[0]]])
    facetFlags3d.append(boundaryTags['screen'])
    facetHoles3d.append([])
    print vertices3d,vertexFlags3d,facets3d,facetFlags3d,facetHoles3d
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices3d,
                                                 vertexFlags=vertexFlags3d,
                                                 facets=facets3d,
                                                 facetFlags=facetFlags3d)
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain = collector3d(fileprefix="collector3dDomainTest",
                        length=1.0,
                        radius=0.1,
                        center=(0.5,0.5,0.5),
                        n_points_on_screen=2*5-2,
                        cross_section=circular_cross_section,
                        thetaOffset=0.0)
    domain.writeAsymptote("collector3dDomainTest")
    domain.writePoly("collector3dDomainTest")
    os.system("asy -V collector3dDomainTest")

