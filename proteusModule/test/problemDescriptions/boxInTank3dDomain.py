import math
from pyadh import Domain

def boxInTank3d(L=[1.0,1.0,1.0],
                box_xy=[0.25,0.25],
                box_L=[0.5,0.5,0.5]):
    boundaries=['left','right','bottom','top','front','back','obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    holes = [[0.5*box_L[0]+box_xy[0],0.5*box_L[1]+box_xy[1],0.5*box_L[2]]]
    vertices=[[0.0,0.0,0.0],#0
              [L[0],0.0,0.0],#1
              [L[0],L[1],0.0],#2
              [0.0,L[1],0.0],#3
              [0.0,0.0,L[2]],#4
              [L[0],0.0,L[2]],#5
              [L[0],L[1],L[2]],#6
              [0.0,L[1],L[2]],#7
              [box_xy[0],box_xy[1],0.0],#8
              [box_xy[0]+box_L[0],box_xy[1],0.0],#9
              [box_xy[0]+box_L[0],box_xy[1]+box_L[1],0.0],#10
              [box_xy[0],box_xy[1]+box_L[1],0.0],#11
              [box_xy[0],box_xy[1],box_L[2]],#12
              [box_xy[0]+box_L[0],box_xy[1],box_L[2]],#13
              [box_xy[0]+box_L[0],box_xy[1]+box_L[1],box_L[2]],#14
              [box_xy[0],box_xy[1]+box_L[1],box_L[2]]]#15
    vertexFlags=[boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left'],
                 boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle'],
                 boundaryTags['obstacle']]
    facets=[[[0,1,2,3],[8,9,10,11]],
            [[0,1,5,4]],
            [[1,2,6,5]],
            [[2,3,7,6]],
            [[3,0,4,7]],
            [[4,5,6,7]],
            [[8,9,13,12]],
            [[9,10,14,13]],
            [[10,11,15,14]],
            [[11,8,12,15]],
            [[12,13,14,15]]]
    facetFlags=[boundaryTags['bottom'],
                boundaryTags['front'],
                boundaryTags['right'],
                boundaryTags['back'],
                boundaryTags['left'],
                boundaryTags['top'],
                boundaryTags['obstacle'],
                boundaryTags['obstacle'],
                boundaryTags['obstacle'],
                boundaryTags['obstacle'],
                boundaryTags['obstacle']]
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
    domain = boxInTank3d()
    domain.writeAsymptote("boxInTank3d")
    domain.writePoly("boxInTank3d")
    domain.writePLY("boxInTank3d")
    os.system("asy -V boxInTank3d")

