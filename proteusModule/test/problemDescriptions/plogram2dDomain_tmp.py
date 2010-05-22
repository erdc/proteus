from pyadh import Domain

def plogram2D(Lx,Ly,slope):
    """
    generate 2D parallelogram domain with boundary flags
    """
    boundaryFlags={'left':1,'back':2,'right':3,'front':4, 'obstacle':5}
    vertices=[]
    vertexFlags=[]
    segments=[]
    segmentFlags=[]

    vertices.append([0.0,0.0])
    vertices.append([0.0, Ly])
    vertices.append([Lx, slope*Lx+Ly])
    vertices.append([Lx,slope*Lx])

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
    domain = plogram2D(10.0, 1.0,-.1)
    domain.writeAsymptote("plogram2D")
    domain.writePoly("plogram2D")
    os.system("asy -V plogram2D")
