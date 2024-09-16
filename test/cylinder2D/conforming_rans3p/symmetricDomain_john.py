import math
from proteus import Domain
import os
import copy

def symmetric2D(box=(1.0,0.41),
                L= 0.15,
                H = 0.15,
                r = 0.05,
                C = (0.15,0.15),
                DX = 0.01,
                refinement_length=0.5,
                DX_coarse = 0.02):
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
    top_right=len(vertices0)-1
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
    nPoints_cyl = max([nPoints_cyl,10])
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
    #section_name = "top_right"
    section_name=os.path.dirname(os.path.abspath(__file__))+"/../conforming_rans2p/"+"top_right"
    #domain0.writePoly(section_name)
    domain0.polyfile=section_name
    genMesh=False
    triangleOptions = "pAq30.0Dena%f" % (.5*DX**2)
    #os.system("triangle " + section_name + ".poly -" + triangleOptions)
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
    groups_hold = copy.copy(groups)
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
        vertices2.append((v[0]+C[0], v[1]+C[1]))

    for i,v in enumerate(vertexFlags):
        if v == boundaryFlags['back']:
            vertexFlags[i] = 0
        elif v == boundaryFlags['right']:
            if vertices2[i][1]<1.0e-8:
                vertexFlags[i]=boundaryFlags['front']
            else:
                vertexFlags[i] = 0
    for i,s in enumerate(segmentFlags):
        if s == boundaryFlags['back']:
            segmentFlags[i] = 0
        elif s == boundaryFlags['right']:
            segmentFlags[i] = 0

    v =len(vertices2)
    vertices2.append((0.0, box[1]))
    vertexFlags.append(boundaryFlags['left'])
    vertices2.append((refinement_length, box[1]))
    vertexFlags.append(boundaryFlags['back'])
    vertices2.append((box[0], box[1]))
    vertexFlags.append(boundaryFlags['right'])
    vertices2.append((box[0], 0.0))
    vertexFlags.append(boundaryFlags['right'])
    vertices2.append((refinement_length, 0.0))
    vertexFlags.append(boundaryFlags['front'])

    segments.append([groups[top_right][1], v])
    segmentFlags.append(boundaryFlags['left'])
    segments.append([v, v+1])
    segmentFlags.append(boundaryFlags['back'])
    segments.append([v+1, v+2])
    segmentFlags.append(boundaryFlags['back'])
    segments.append([v+2, v+3])
    segmentFlags.append(boundaryFlags['right'])
    segments.append([v+3, v+4])
    segmentFlags.append(boundaryFlags['front'])
    segments.append([v+4,groups_hold[top_right][1]])
    segmentFlags.append(boundaryFlags['front'])

    segments.append([v+1,v+4])
    segmentFlags.append(0)

    regions0 = [[vertices2[v][0]+1.0e-8, vertices2[v][1]-1.0e-8],
                [vertices2[v+2][0]-1.0e-8, vertices2[v+2][1]-1.0e-8]]
    regionFlags0=[0,0]
    regionConstraints = [0.5*DX**2, 0.5*DX_coarse**2]
    
    
            
    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices2,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  holes=[C],
                                                  regions=regions0,
                                                  regionFlags=regionFlags0,
                                                  regionConstraints=regionConstraints)
    domain.boundaryFlags = boundaryFlags
    return domain
    #domain1.writePoly("test")


 
  
    
    
    


if __name__=='__main__':
    import os
    domain = symmetric2D()
    domain.writePoly("symmetric2D")
    
 
        
        
