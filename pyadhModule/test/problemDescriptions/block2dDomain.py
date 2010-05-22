import math
from pyadh import Domain

def block2D(block_size,
            points_on_grain=50,
            points_on_boundary=50):
    """
    generates a rectagular domain of a given size    """
    #n_domain_vertices = 4
    #domain_vertices =[(0.0,0.0),(0.0,1.0),(1.0,1.0),(1.0,0.0)]
    n_domain_vertices = 4 + 4*points_on_boundary
    DX=block_size[0]/float(points_on_boundary+1)
    DY=block_size[1]/float(points_on_boundary+1)
    domain_vertices=[]
    for k in range(points_on_boundary+2):
        domain_vertices.append((0.0,k*DY))
    for k in range(1,points_on_boundary+2):
        domain_vertices.append((k*DX,block_size[1]))
    for k in range(points_on_boundary,-1,-1):
        domain_vertices.append((block_size[0],k*DY))
    for k in range(points_on_boundary,0,-1):
        domain_vertices.append((k*DX,0.0))

    nvertices = len(domain_vertices)
    boundaries = ['left', 'right', 'front', 'back', 'obstacle']
    #based on vertex order and way initial segment list goes round
    boundaryFlags={'left':1,'back':2,'right':3,'front':4, 'obstacle':5}
    vertices=[]
    vertexFlags=[]
    segments=[]
    segmentFlags=[]
    #write vertices
    
    for v,p in enumerate(domain_vertices):#numbering is base 1 for triangle
        if (p[0] == 0.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['left'])
        elif (p[0] == block_size[0]):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['right'])
        elif (p[1] == 0.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['front'])
        elif (p[1] == block_size[1]):
            vertexFlags.append(boundaryFlags['back'])
            vertices.append([p[0],p[1]])
        else:
            exit
    #write segments
    #left is X_minus, right is X_plus, front is Y_minus, back is Y_plus
    segments=[]
    segmentFlags=[]
    for sN in range(len(domain_vertices)-1):
        segments.append([sN,sN+1])
        if (domain_vertices[sN][0] == 0.0 and domain_vertices[sN+1][0] ==  0.0):
            segmentFlags.append(boundaryFlags['left'])
        elif (domain_vertices[sN][0] == block_size[0] and domain_vertices[sN+1][0] ==  block_size[0]):
            segmentFlags.append(boundaryFlags['right'])
        elif (domain_vertices[sN][1] == 0.0 and domain_vertices[sN+1][1] ==  0.0):
            segmentFlags.append(boundaryFlags['front'])
        elif (domain_vertices[sN][1] == block_size[1] and domain_vertices[sN+1][1] ==  block_size[1]):
            segmentFlags.append(boundaryFlags['back'])
        else:
            exit
    segments.append([len(domain_vertices)-1,0])
    if (domain_vertices[segments[-1][0]][0] == 0.0 and domain_vertices[segments[-1][1]][0] ==  0.0):
        segmentFlags.append(boundaryFlags['left'])
    if (domain_vertices[segments[-1][0]][0] == block_size[0] and domain_vertices[segments[-1][1]][0] ==  block_size[0]):
        segmentFlags.append(boundaryFlags['right'])
    if (domain_vertices[segments[-1][0]][1] == 0.0 and domain_vertices[segments[-1][1]][1] ==  0.0):
        segmentFlags.append(boundaryFlags['front'])
    if (domain_vertices[segments[-1][0]][1] == block_size[1] and domain_vertices[segments[-1][1]][1] ==  block_size[1]):
        segmentFlags.append(boundaryFlags['back'])
    #end exterior boundary segments

    holes=[]

    
    regions=[[vertices[0][0]+1.0e-8,
              vertices[0][1]+1.0e-8]]
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
    domain.boundaryFlags = boundaryFlags

    return domain

if __name__=='__main__':
    import os
    domain = block2D(block_size= [2.0, 0.5],points_on_grain=50, points_on_boundary=5)
    domain.writeAsymptote("block2D")
    domain.writePoly("block2D")
    os.system("asy -V block2D")
