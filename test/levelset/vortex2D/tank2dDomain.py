import math
from proteus import Domain

def tank2d(L=[1.0,1.0,1.0],fileprefix=None):
    boundaries=['left','right','bottom','top','front','back','obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    vertices=[[0.0,0.0],#0
              [L[0],0.0],#1
              [L[0],L[1]],#2
              [0.0,L[1]]]#3
    vertexFlags=[boundaryTags['left'],
                 boundaryTags['right'],
                 boundaryTags['right'],
                 boundaryTags['left']]
    segments=[[0,1],
              [1,2],
              [2,3],
              [3,0]]
    segmentFlags=[boundaryTags['front'],
                  boundaryTags['right'],
                  boundaryTags['back'],
                  boundaryTags['left']]
    regions=[[0.5*L[0],0.5*L[1]]]
    regionFlags=[1.0]
    domain = Domain.PlanarStraightLineGraphDomain(fileprefix=fileprefix,
                                                  vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  regions=regions,
                                                  regionFlags=regionFlags)
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain = tank2d()
    domain.writePoly("tank2d")
    domain.writePLY("tank2d")
    os.system("asy -V tank2d")
