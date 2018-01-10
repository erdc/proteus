import math
from proteus import Domain

def cylinder2D(L=2.5,
               H=0.41,
               C=(0.2,0.2),
               r=0.05,
               points_on_grain=10):
    """
    makes the cylinder2d  domain.
    """
    
    domain_vertices = [(0.0,0.0), (0.0,H),  (L,H),(L,0.0)]
    vertices = []
    vertexFlags = []
    segments = []
    segmentFlags = []
    grain_centers = [C] #[(0.2,0.2)]
    radius = r
    boundaries = ['left', 'right', 'front', 'back', 'obstacle']
    boundaryFlags=dict([(key,i+1) for i,key in enumerate(boundaries)])

    

    for v,p in enumerate(domain_vertices):
        if (p[0] == 0.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['left'])
        elif (p[0] == L):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['right'])
        elif (p[1] == 0.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['front'])
        elif (p[1] == H):
            vertexFlags.append(boundaryFlags['back'])
            vertices.append([p[0],p[1]])
        else:
            exit
    for sN in range(len(domain_vertices)-1):
        segments.append([sN,sN+1])
        if (domain_vertices[sN][0] == 0.0 and domain_vertices[sN+1][0] ==  0.0):
            segmentFlags.append(boundaryFlags['left'])
        elif (domain_vertices[sN][0] == L and domain_vertices[sN+1][0] ==  L):            segmentFlags.append(boundaryFlags['right'])
        elif (domain_vertices[sN][1] == 0.0 and domain_vertices[sN+1][1] ==  0.0):
            segmentFlags.append(boundaryFlags['front'])
        elif (domain_vertices[sN][1] == H and domain_vertices[sN+1][1] ==  H):
            segmentFlags.append(boundaryFlags['back'])
        else:
            exit
    segments.append([len(domain_vertices)-1,0])
    if (domain_vertices[segments[-1][0]][0] == 0.0 and domain_vertices[segments[-1][1]][0] ==  0.0):
        segmentFlags.append(boundaryFlags['left'])
    if (domain_vertices[segments[-1][0]][0] == L and domain_vertices[segments[-1][1]][0] ==  L):
        segmentFlags.append(boundaryFlags['right'])
    if (domain_vertices[segments[-1][0]][1] == 0.0 and domain_vertices[segments[-1][1]][1] ==  0.0):
        segmentFlags.append(boundaryFlags['front'])
    if (domain_vertices[segments[-1][0]][1] == H and domain_vertices[segments[-1][1]][1] ==  H):
        segmentFlags.append(boundaryFlags['back'])

    vStart = len(domain_vertices)
    sStart = len(segments)
    for g,c in enumerate(grain_centers):
        for gb in range(points_on_grain):
            vertices.append([c[0]+radius*math.sin(float(gb)/float(points_on_grain)*2.0*math.pi),c[1]+radius*math.cos(float(gb)/float(points_on_grain)*2.0*math.pi)])
            vertexFlags.append(boundaryFlags['obstacle'])
    for rb in range(len(grain_centers)):
        for gb in range(points_on_grain-1):
            segments.append([sStart+points_on_grain*rb+gb,sStart+points_on_grain*rb+gb+1])
            segmentFlags.append(boundaryFlags['obstacle'])
        segments.append([sStart+points_on_grain*rb+points_on_grain-1,sStart+points_on_grain*rb])
        segmentFlags.append(boundaryFlags['obstacle'])

   
    
    regions=[[vertices[0][0]+1.0e-8,
              vertices[0][1]+1.0e-8]]
    regionFlags=[1]
    


    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  holes=grain_centers,
                                                  regions=regions,
                                                  regionFlags=regionFlags)
    #go ahead and add a boundary tags member 
    domain.boundaryFlags = boundaryFlags

    return domain

if __name__=='__main__':
    import os
    domain = cylinder2D()
    #domain.writeAsymptote("pome2D")
    domain.writePoly("cylinder2D")
    domain.writePLY("cylinder2D")
