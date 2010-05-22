import math
from pyadh import Domain

def wavetank3d(water_sponge_length=1.0,
            piston_length=1.0,
            bottom_length=1.0,
            beach_slope=2.0,
            beach_length_x=1.0,
            land_length=1.0,
            land_sponge_length=1.0,
            back_height=1.0,
            width=1.0,
            attack_slope=1.0):
    boundaries=['piston',
                'walls',
                'top']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    vertices_front=[(0.0,
                     0.0,
                     0.0,
                     boundaryTags['walls']),
                    (water_sponge_length,
                     0.0,
                     0.0,
                     boundaryTags['walls']),
                    (water_sponge_length+piston_length,
                     0.0,
                     0.0,
                     boundaryTags['walls']),
                    (water_sponge_length+piston_length+bottom_length,
                     0.0,
                     0.0,
                     boundaryTags['walls']),
                    (water_sponge_length+piston_length+bottom_length+beach_length_x,
                     0.0,
                     beach_slope*beach_length_x,
                     boundaryTags['walls']),
                    (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length,
                     0.0,
                     beach_slope*beach_length_x,
                     boundaryTags['walls']),
                    (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+land_sponge_length,
                     0.0,
                     beach_slope*beach_length_x,
                     boundaryTags['walls']),
                    (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+land_sponge_length,
                     0.0,
                     beach_slope*beach_length_x+back_height,
                     boundaryTags['walls']),
                                        (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length,
                                         0.0,
                                         beach_slope*beach_length_x+back_height,
                                         boundaryTags['walls']),
                    (water_sponge_length,
                     0.0,
                     beach_slope*beach_length_x+back_height,
                     boundaryTags['walls']),
                    (0.0,
                     0.0,
                     beach_slope*beach_length_x+back_height,
                     boundaryTags['walls'])]
    vertices_back=[(0.0,
                    width,
                    0.0,
                    boundaryTags['walls']),
                   (water_sponge_length,
                    width,
                    0.0,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length,
                    width,
                    0.0,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length+bottom_length+attack_slope,
                    width,
                    0.0,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length+bottom_length+beach_length_x+attack_slope,
                    width,
                    beach_slope*beach_length_x,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+attack_slope,
                    width,
                    beach_slope*beach_length_x,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+land_sponge_length+attack_slope,
                    width,
                    beach_slope*beach_length_x,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+land_sponge_length+attack_slope,
                    width,
                    beach_slope*beach_length_x+back_height,
                    boundaryTags['walls']),
                   (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+attack_slope,
                    width,
                    beach_slope*beach_length_x+back_height,
                    boundaryTags['walls']),
                   (water_sponge_length,
                    width,
                    beach_slope*beach_length_x+back_height,
                    boundaryTags['walls']),
                   (0.0,
                    width,
                    beach_slope*beach_length_x+back_height,
                    boundaryTags['walls'])]
    nVertices_front=len(vertices_front)
    facets=[[[1,2,10,11],boundaryTags['walls']],
            [[2,3,4,5,6,9,10],boundaryTags['walls']],
            [[6,7,8,9],boundaryTags['walls']],
            [[12,13,21,22],boundaryTags['walls']],
            [[13,14,15,16,17,20,21],boundaryTags['walls']],
            [[17,18,19,20],boundaryTags['walls']]]
    for vN in range(nVertices_front):
        if vN >= 7 and vN <= 9:
            tag = boundaryTags['top']
        else:
            tag = boundaryTags['walls']
        facets.append([[1+vN,
                        1+(vN+1)%nVertices_front,
                        1+(vN+1)%nVertices_front+nVertices_front,
                        1+vN+nVertices_front],
                       tag])
    facets.append([[2,nVertices_front-1,nVertices_front+nVertices_front-1,nVertices_front+2],0])
    facets.append([[1+nVertices_front-6,1+nVertices_front-3,1+nVertices_front+nVertices_front-3,1+nVertices_front+nVertices_front-6],0])
    #label piston facet
    facets[7][1]=boundaryTags['piston']
    #rewrite to new format
    newVertices=[]
    newVertexFlags=[]
    newFacets=[]
    newFacetFlags=[]
    for v in vertices_front+vertices_back:
        newVertices.append([v[0],v[1],v[2]])        
        newVertexFlags.append(v[3])
    for f in facets:
        newFacets.append([[vN-1 for vN in f[0]]])        
        newFacetFlags.append(f[1])
    print newVertices
    print newVertexFlags
    print newFacets
    print newFacetFlags
    domain = Domain.PiecewiseLinearComplexDomain(vertices=newVertices,
                                                 vertexFlags=newVertexFlags,
                                                 facets=newFacets,
                                                 facetFlags=newFacetFlags,
                                                 regions=[[newVertices[0][0]+1.0e-8,
                                                          newVertices[0][1]+1.0e-8,
                                                           newVertices[0][2]+1.0e-8]],
                                                 regionFlags=[1])
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain =  wavetank3d(water_sponge_length=1.0,
                         piston_length=1.0,
                         bottom_length=1.0,
                         beach_slope=2.0,
                         beach_length_x=1.0,
                         land_length=1.0,
                         land_sponge_length=1.0,
                         back_height=1.0,
                         width=1.0,
                         attack_slope=1.0)
    domain.writeAsymptote("wavetank3d")
    domain.writePoly("wavetank3d")
    domain.writePLY("wavetank3d")
    os.system("asy -V wavetank3d")

