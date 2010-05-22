

def genPoly(fileprefix,
            water_sponge_length=1.0,
            piston_length=1.0,
            bottom_length=1.0,
            beach_slope=2.0,
            beach_length_x=1.0,
            land_length=1.0,
            land_sponge_length=1.0,
            back_height=1.0):
    boundaries=['piston',
                'walls',
                'top']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    vertices=[(0.0,0.0,boundaryTags['walls']),
              (water_sponge_length,0.0,boundaryTags['walls']),
              (water_sponge_length+piston_length,0.0,boundaryTags['walls']),
              (water_sponge_length+piston_length+bottom_length,0.0,boundaryTags['walls']),
              (water_sponge_length+piston_length+bottom_length+beach_length_x,
               beach_slope*beach_length_x,
               boundaryTags['walls']),
              (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length,
               beach_slope*beach_length_x,
               boundaryTags['walls']),
              (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+land_sponge_length,
               beach_slope*beach_length_x,
               boundaryTags['walls']),
              (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length+land_sponge_length,
               beach_slope*beach_length_x+back_height,
               boundaryTags['walls']),
              (water_sponge_length+piston_length+bottom_length+beach_length_x+land_length,
               beach_slope*beach_length_x+back_height,
               boundaryTags['walls']),
              (water_sponge_length,
               beach_slope*beach_length_x+back_height,
               boundaryTags['walls']),
              (0.0,
               beach_slope*beach_length_x+back_height,
               boundaryTags['walls'])]
    segments=[(1,2,boundaryTags['walls']),
              (2,3,boundaryTags['piston']),
              (3,4,boundaryTags['walls']),
              (4,5,boundaryTags['walls']),
              (5,6,boundaryTags['walls']),
              (6,7,boundaryTags['walls']),
              (7,8,boundaryTags['walls']),
              (8,9,boundaryTags['top']),
              (9,10,boundaryTags['top']),
              (10,11,boundaryTags['top']),
              (11,1,boundaryTags['walls']),
              (2,10,0),
              (6,9,0)]
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
    poly.write("#vertices \n")
    for vN,v in enumerate(vertices):
        poly.write('%d %12.5e %12.5e %d\n' % (vN+1,v[0],v[1],v[2]))
    poly.write("#segments \n")
    poly.write('%d %d \n' % (len(segments),1))
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0],s[1],s[2]))
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (2,))
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          0.5*(vertices[0][0]+vertices[1][0]),
                                          0.5*(vertices[0][1]+vertices[-1][1]),
                                          1))
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          0.5*(vertices[-3][0]+vertices[-4][0]),
                                          0.5*(vertices[-4][1]+vertices[-5][1]),
                                          2))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('wavetank2d_test')
    os.system('triangle -Aq30den wavetank2d_test.poly')
    os.system('showme wavetank2d_test.1.ele')
