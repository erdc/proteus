import math
from pyadh import Domain
import random

def boxesInTank3d(L=[1.0,1.0,1.0],
                nBoxes_xy=[10,10],
                boxScale=[0.5,0.5,0.35]):
    rg = random.Random()
    boundaries=['left','right','bottom','top','front','back','obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    vertices=[[0.0,0.0,0.0],#0
              [L[0],0.0,0.0],#1
              [L[0],L[1],0.0],#2
              [0.0,L[1],0.0],#3
              [0.0,0.0,L[2]],#4
              [L[0],0.0,L[2]],#5
              [L[0],L[1],L[2]],#6
              [0.0,L[1],L[2]]]#7
#               [box_xy[0],box_xy[1],0.0],#8
#               [box_xy[0]+box_L[0],box_xy[1],0.0],#9
#               [box_xy[0]+box_L[0],box_xy[1]+box_L[1],0.0],#10
#               [box_xy[0],box_xy[1]+box_L[1],0.0],#11
#               [box_xy[0],box_xy[1],box_L[2]],#12
#               [box_xy[0]+box_L[0],box_xy[1],box_L[2]],#13
#               [box_xy[0]+box_L[0],box_xy[1]+box_L[1],box_L[2]],#14
#               [box_xy[0],box_xy[1]+box_L[1],box_L[2]]]#15
    vertexFlags=[boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left'],
                 boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left']]
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle'],
#                  boundaryTags['obstacle']]
    facets=[[[0,1,2,3]],
            [[0,1,5,4]],
            [[1,2,6,5]],
            [[2,3,7,6]],
            [[3,0,4,7]],
            [[4,5,6,7]]]
#             [[8,9,13,12]],
#             [[9,10,14,13]],
#             [[10,11,15,14]],
#             [[11,8,12,15]],
#             [[12,13,14,15]]]
    facetFlags=[boundaryTags['bottom'],
                boundaryTags['front'],
                boundaryTags['right'],
                boundaryTags['back'],
                boundaryTags['left'],
                boundaryTags['top']]
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle'],
#                 boundaryTags['obstacle']]
    holes=[]
    dx = L[0]/float(nBoxes_xy[0])
    dy = L[1]/float(nBoxes_xy[1])
    hbdx = 0.5*boxScale[0]*dx
    hbdy = 0.5*boxScale[1]*dy
    for i in range(nBoxes_xy[0]):
        for j in range(nBoxes_xy[1]):
            center = [i*dx + 0.5*dx,j*dy+0.5*dy]
            holes.append(center+[0.5*boxScale[2]*L[2]])
            hz = rg.gauss(1.0,0.1)*boxScale[2]*L[2]
            n=len(vertices)
            vertices.append([center[0]-hbdx,center[1]-hbdy,0.0])
            vertices.append([center[0]-hbdx,center[1]+hbdy,0.0])
            vertices.append([center[0]+hbdx,center[1]+hbdy,0.0])
            vertices.append([center[0]+hbdx,center[1]-hbdy,0.0])
            vertices.append([center[0]-hbdx,center[1]-hbdy,hz])
            vertices.append([center[0]-hbdx,center[1]+hbdy,hz])
            vertices.append([center[0]+hbdx,center[1]+hbdy,hz])
            vertices.append([center[0]+hbdx,center[1]-hbdy,hz])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            vertexFlags.append(boundaryTags['obstacle'])
            facets[0].append([n+0,n+1,n+2,n+3])#bottom
            facets.append([[n+0,n+1,n+5,n+4]])#left
            facets.append([[n+1,n+2,n+6,n+5]])#back
            facets.append([[n+2,n+3,n+7,n+6]])#right
            facets.append([[n+3,n+0,n+4,n+7]])#front
            facets.append([[n+4,n+5,n+6,n+7]])#top
            facetFlags.append(boundaryTags['obstacle'])
            facetFlags.append(boundaryTags['obstacle'])
            facetFlags.append(boundaryTags['obstacle'])
            facetFlags.append(boundaryTags['obstacle'])
            facetFlags.append(boundaryTags['obstacle'])
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 vertexFlags=vertexFlags,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 holes=holes)
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain = boxesInTank3d()
    domain.writeAsymptote("boxesInTank3d")
    domain.writePoly("boxesInTank3d")
    domain.writePLY("boxesInTank3d")
    os.system("asy -V boxesInTank3d")

