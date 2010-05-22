import math
from pyadh import Domain

def vesselPoly3d(fileprefix,
                 vesselprefix,
                 outerBoxFactor=[0.25,0.25,0.25],
                 offsetFactor=[0.6,0.0,0.0]):
    vesselDomain = Domain.PiecewiseLinearComplexDomain(units="m")
    vesselDomain.readPoly(vesselprefix)
    for vN,v in enumerate(vesselDomain.vertices):
        vesselDomain.vertices[vN]=[v[0]*0.3048,v[1]*0.3048,v[2]*0.3048,]
    vesselDomain.getBoundingBox()
    boundaries=['left',
                'right',
                'front',
                'back',
                'bottom',
                'top',
                'obstacle']
    L = [l + f*vesselDomain.L[1] for l,f in zip(vesselDomain.L,outerBoxFactor)]
    x = [xVessel-0.5*f*vesselDomain.L[1] + of*0.5*l for xVessel,f,of,l in zip(vesselDomain.x,outerBoxFactor,offsetFactor,L)]
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    nVertices_vessel = len(vesselDomain.vertices)
    vesselDomain.vertexFlags=[boundaryTags['obstacle'] for v in vesselDomain.vertices]
    vesselDomain.facetFlags=[boundaryTags['obstacle'] for f in vesselDomain.facets]
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [[x[0], x[1],x[2]],#0
                [x[0] + L[0],x[1],x[2]],#1
                [x[0] + L[0],x[1],x[2]+L[2]],#2
                [x[0], x[1],x[2]+L[2]],#3
                [x[0], x[1]+L[1],x[2]],#4
                [x[0] + L[0],x[1]+L[1],x[2]],#5
                [x[0] + L[0],x[1]+L[1],x[2]+L[2]],#6
                [x[0], x[1]+L[1],x[2]+L[2]]]#7
    vertexFlags = [boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top'],
                   boundaryTags['bottom'],
                   boundaryTags['bottom'],
                   boundaryTags['top'],
                   boundaryTags['top']]
    facets =[[[nVertices_vessel+0,nVertices_vessel+1,nVertices_vessel+2,nVertices_vessel+3]],#front
             [[nVertices_vessel+4,nVertices_vessel+5,nVertices_vessel+6,nVertices_vessel+7]],#back
             [[nVertices_vessel+0,nVertices_vessel+1,nVertices_vessel+5,nVertices_vessel+4]],#bottom
             [[nVertices_vessel+2,nVertices_vessel+3,nVertices_vessel+7,nVertices_vessel+6]],#top
             [[nVertices_vessel+0,nVertices_vessel+3,nVertices_vessel+7,nVertices_vessel+4]],#left
             [[nVertices_vessel+1,nVertices_vessel+2,nVertices_vessel+6,nVertices_vessel+5]]]#right
    facetFlags=[boundaryTags['front'],
                boundaryTags['back'],
                boundaryTags['bottom'],
                boundaryTags['top'],
                boundaryTags['left'],
                boundaryTags['right']]
    print L
    print x
    holes=[[vesselDomain.x[0]+0.5*vesselDomain.L[0],vesselDomain.x[1]+0.5*vesselDomain.L[1],vesselDomain.x[2]+0.5*vesselDomain.L[2]]]
    def newFlag(f):
        if f == boundaryTags['obstacle']:
           return f
        else:
           return 0
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vesselDomain.vertices + vertices,
                                                 vertexFlags=[newFlag(vf) for vf in vesselDomain.vertexFlags]+vertexFlags,
                                                 facets=vesselDomain.facets+facets,
                                                 facetFlags=[newFlag(vf) for ff in vesselDomain.facetFlags]+facetFlags,
                                                 holes=holes)
#    domain = Domain.PiecewiseLinearComplexDomain(vertices=vesselDomain.vertices + vertices,
#                                                 vertexFlags=vesselDomain.vertexFlags+vertexFlags,
#                                                 facets=vesselDomain.facets+facets,
#                                                 facetFlags=vesselDomain.facetFlags+facetFlags,
#                                                 holes=holes)
#     domain = Domain.PiecewiseLinearComplexDomain(vertices=vesselDomain.vertices,
#                                                  vertexFlags=vesselDomain.vertexFlags,
#                                                  facets=vesselDomain.facets,
#                                                  facetFlags=vesselDomain.facetFlags,
#                                                  holes=holes)
    domain.boundaryTags = boundaryTags
    domain.hull_length = vesselDomain.L[0]
    domain.hull_beam = vesselDomain.L[1]
    domain.hull_draft = vesselDomain.L[2]
    domain.getBoundingBox()
    return domain

if __name__=='__main__':
    import os
    domain = vesselPoly3d(fileprefix="vesselPoly3dDomainTest",
                          vesselprefix="Container06_3d",
                          outerBoxFactor=[10.0,10.0,0.1])
    domain.writePoly("vesselPoly3dDomainTest")
    domain.writePLY("vesselPoly3dDomainTest")
    os.system("tetgen -KVApq1.25q13fen vesselPoly3dDomainTest")
    #domain.writeAsymptote("wigley3dDomainTest")

