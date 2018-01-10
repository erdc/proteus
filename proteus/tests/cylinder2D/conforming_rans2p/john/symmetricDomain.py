import math
from proteus import Domain
import os
import copy

def symmetric2D(L= 1.0,
                H = 1.0,
                r = 0.5,
                DX = 0.2):
    boundaries = ['left', 'right', 'front', 'back', 'obstacle']
    boundaryFlags=dict([(key,i+1) for i,key in enumerate(boundaries)])
    left_nPoints = int(math.ceil((H-r)/DX))
    top_nPoints = int(math.ceil(L/DX))
    right_nPoints = int(math.ceil(H/DX))
    bottom_nPoints = int(math.ceil((L-r)/DX))
    DX_left = (H-r)/float(left_nPoints)
    DX_top = L/float(top_nPoints)
    DX_right = H/float(right_nPoints)
    DX_bottom = (L-r)/float(bottom_nPoints)

    vertices0 = [(0.0, r)]
    vertexFlags0 = [boundaryFlags['obstacle']]
    segments0 = []
    segmentFlags0 = []
    grain_centers = [(0.0,0.0)] #[(0.2,0.2)]
    radius = r

    # left
    for i in range(1,left_nPoints):
        vertices0.append((0.0, r+i*DX_left))
        vertexFlags0.append(0)
    # top
    for i in range(0,top_nPoints+1):
        vertices0.append((i*DX_top, H))
        vertexFlags0.append(boundaryFlags['back'])
    # right
    for i in range(1, right_nPoints+1):
        vertices0.append((L, H-i*DX_right))
        vertexFlags0.append(boundaryFlags['right'])
    # bottom
    for i in range(1, bottom_nPoints):
        vertices0.append((L-i*DX_bottom, 0.0))
        vertexFlags0.append(0)

    arclength= 0.5*math.pi*r
    nPoints_cyl = int(math.ceil(arclength/DX))
    #DX_cyl = arclength/float(nPoints_cyl)

    # cyl
    for i in range(nPoints_cyl):
        vertices0.append((r*math.cos(float(i)/float(nPoints_cyl)*0.5*math.pi),r*math.sin(float(i)/float(nPoints_cyl)*0.5*math.pi)))
        vertexFlags0.append(boundaryFlags['obstacle'])

    for sN in range(len(vertices0)-1):
        segments0.append([sN,sN+1])
        if (vertices0[sN][0] == 0.0 and vertices0[sN+1][0] ==  0.0):
            segmentFlags0.append(0)
        elif (vertices0[sN][0] == L and vertices0[sN+1][0] ==  L):            
            segmentFlags0.append(boundaryFlags['right'])
        elif (vertices0[sN][1] == 0.0 and vertices0[sN+1][1] ==  0.0):
            segmentFlags0.append(0)
        elif (vertices0[sN][1] == H and vertices0[sN+1][1] ==  H):
            segmentFlags0.append(boundaryFlags['back'])
        else:
            segmentFlags0.append(boundaryFlags['obstacle'])
    
    segments0.append([len(vertices0)-1,0])    
    segmentFlags0.append(boundaryFlags['obstacle'])

    regions0=[[L-1.0e-8,
              H-1.0e-8]]
    regionFlags0=[0]
    domain0 = Domain.PlanarStraightLineGraphDomain(vertices=vertices0,
                                                  vertexFlags=vertexFlags0,
                                                  segments=segments0,
                                                  segmentFlags=segmentFlags0,
                                                  holes=[],
                                                  regions=regions0,
                                                  regionFlags=regionFlags0)
    section_name = "top_right"
    domain0.writePoly(section_name)
    triangleOptions = "pAq30.0Dena%f" % (.5*DX**2)
    os.system("triangle " + section_name + ".poly -" + triangleOptions)
    nodes = open(section_name + ".1.node", 'r')
    edges = open(section_name + ".1.edge", 'r')



    vertices = []
    vertexFlags = []
    segments = []
    segmentFlags = []


    a = nodes.readline().split()
    nodeCount =int(a[0])
    for i in range(1, nodeCount+1):
        b= nodes.readline().split()
        vertices.append((float(b[1]),float(b[2])))
        if int(b[3])==1:
            vertexFlags.append(0)
        else:
            vertexFlags.append(int(b[3]))
    nodes.close()

    a = edges.readline().split()
    nodeCount =int(a[0])
    for i in range(1, nodeCount+1):
        b= edges.readline().split()
        segments.append([int(b[1])-1,int(b[2])-1])
        if int(b[3])==1:
            segmentFlags.append(0)
        else:
            segmentFlags.append(int(b[3]))
    edges.close()

    groups=[]
    vertices_hold = copy.copy(vertices)
    for i,v in enumerate(vertices_hold):
        groups.append([i])
        if v[1]==0.0:
            groups[i].append(i)
        else:
            vertices.append((v[0],-v[1]))
            if vertexFlags[i] == boundaryFlags['back']:
                vertexFlags.append(boundaryFlags['front'])
            else:
                vertexFlags.append(vertexFlags[i])
            groups[i].append(len(vertices)-1)
    segments_hold = copy.copy(segments)

    for i,s in enumerate(segments_hold):
        sN = [groups[s[0]][1], groups[s[1]][1]]
        if sN not in segments:
            segments.append(sN)
            if segmentFlags[i]==boundaryFlags['back']:
                segmentFlags.append(boundaryFlags['front'])
            else:
                segmentFlags.append(segmentFlags[i])
    ###
    groups=[]
    vertices_hold = copy.copy(vertices)
    for i,v in enumerate(vertices_hold):
        groups.append([i])
        if v[0]==0.0:
            groups[i].append(i)
        else:
            vertices.append((-v[0],v[1]))
            if vertexFlags[i] == boundaryFlags['right']:
                vertexFlags.append(boundaryFlags['left'])
            else:
                vertexFlags.append(vertexFlags[i])
            groups[i].append(len(vertices)-1)
    segments_hold = copy.copy(segments)

    for i,s in enumerate(segments_hold):
        sN = [groups[s[0]][1], groups[s[1]][1]]
        if sN not in segments:
            segments.append(sN)
            if segmentFlags[i]==boundaryFlags['right']:
                segmentFlags.append(boundaryFlags['left'])
            else:
                segmentFlags.append(segmentFlags[i])
    ###

    vertices2=[]
    for v in vertices:
        vertices2.append((v[0]+L, v[1]+H))

    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices2,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  holes=[(L,H)],
                                                  regions=regions0,
                                                  regionFlags=regionFlags0)
    domain.boundaryFlags = boundaryFlags
    return domain
    #domain1.writePoly("test")


 
  
    
    
    


if __name__=='__main__':
    import os
    domain = symmetric2D()
    domain.writePoly("symmetric2D")
    
 
        
        
