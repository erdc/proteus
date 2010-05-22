

def genPoly(fileprefix,
            water_sponge_length=1.0,
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
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (2*nVertices_front,3,0,1))
    poly.write("#vertices \n")
    for vN,v in enumerate(vertices_front+vertices_back):
        poly.write('%d %12.5e %12.5e %12.5e %d\n' % (vN+1,v[0],v[1],v[2],v[3]))
    poly.write("#facets \n")
    poly.write('%d %d \n' % (len(facets),1))
    for fN,f in enumerate(facets):
        poly.write('1 0 %d\n' % (f[1],))
        poly.write('%d ' % (len(f[0]),))
        for vN in f[0]:
            poly.write('%d ' % (vN,))
        poly.write('\n')
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (3,))
    poly.write('%d %12.5e %12.5e %12.5e 1\n' % (1,
                                               0.5*(vertices_front[0][0]+vertices_front[1][0]),
                                               0.5*(vertices_front[0][1]+vertices_back[0][1]),
                                               0.5*(vertices_front[0][2]+vertices_front[-1][2])))
    poly.write('%d %12.5e %12.5e %12.5e 2\n' % (2,
                                                0.5*(vertices_front[5][0]+vertices_front[6][0]),
                                                0.5*(vertices_front[5][1]+vertices_back[5][1]),
                                                0.5*(vertices_front[5][2]+vertices_front[8][2])))
    poly.write('%d %12.5e %12.5e %12.5e 3\n' % (3,
                                                0.5*(vertices_front[2][0]+vertices_front[3][0]),
                                                0.5*(vertices_front[2][1]+vertices_back[2][1]),
                                                0.5*(vertices_front[2][2]+vertices_front[9][2])))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('wavetank3d_test')
    os.system('tetgen -Aq1.25en wavetank3d_test.poly')
    os.system('tetview wavetank3d_test.1.ele')
