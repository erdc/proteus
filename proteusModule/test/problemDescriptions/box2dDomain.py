from pyadh import Domain

def box2D(Lx,Ly):
    """
    generate 2D rectangular domain with boundary flags
    """
    boundaryFlags={'left':1,'back':2,'right':3,'front':4, 'obstacle':5}
    vertices=[]
    vertexFlags=[]
    segments=[]
    segmentFlags=[]

    vertices.append([0.0,0.0])
    vertices.append([0.0, Ly])
    vertices.append([Lx, Ly])
    vertices.append([Lx, 0.0])

    vertexFlags.append(boundaryFlags['front'])
    vertexFlags.append(boundaryFlags['back'])
    vertexFlags.append(boundaryFlags['back'])
    vertexFlags.append(boundaryFlags['front'])

    for i in range(3):
        segments.append([i, i+1])

    segments.append([3,0])
    segmentFlags.append(boundaryFlags['left'])
    segmentFlags.append(boundaryFlags['back'])
    segmentFlags.append(boundaryFlags['right'])
    segmentFlags.append(boundaryFlags['front'])

    regions=[[vertices[0][0]+1.0e-8,
              vertices[0][1]+1.0e-8]]
    regionFlags=[1]
    holes=[]

    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  holes=holes,
                                                  regions=regions,
                                                  regionFlags=regionFlags)
    domain.boundaryFlags = boundaryFlags
    
    return domain

if __name__=='__main__':
    import os
    domain = box2D(10.0, 1.0)
    domain.writeAsymptote("box2D")
    domain.writePoly("box2D")
    os.system("asy -V box2D")
