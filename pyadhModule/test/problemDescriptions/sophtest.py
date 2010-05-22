import math

#This is the channel from the side
#              top
#             ____________
#            |  _         |        
# upstream   | | |        |  downstream
#            |  -         |
#             ____________
#              bottom
#^
#y x>

def circular_cross_section(center,radius,thetat):

    x = radius*math.cos(theta)+center
    y = radius*math.sin(theta)+center
    return (x,y)

def genPoly(fileprefix,
            height=4.0,
            length_to_height_ration=8.0,
            n_points_on_obstacle=15.0,
            cross_section=circular_cross_section,
            center=(2.5,1.5)):
    boundaries=['upstream',
                'downstream',
                'bottom',
                'top',
                'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    length = height*length_to_height_ration
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices= [('upstream_bottom', (0.0, 0.0),boundaryTags['bottom']),
               ('downstream_bottom',(length,0.0),boundaryTags['bottom']),
               ('downstream_top',(length,height),boundaryTags['top']),
               ('upstream_top',(0.0,height),boundaryTags['top'])]
    nv = len(vertices)
    segments = [(1,2,boundaryTags['bottom']),
                (2,3,boundaryTags['downstream']),
                (3,4,boundaryTags['top']),
                (4,1,boundaryTags['upstream'])]
    #add vertices for square obstacle
    vertices+=[('lower_lefthand_corner',(2.0,1.0),boundaryTags['obstacle']),
               ('lower_righthand_corner',(3.0,1.0),boundaryTags['obstacle']),
               ('upper_righthand_corner',(3.0,2.0),boundaryTags['obstacle']),
               ('upper_lefthand_corner',(2.0,2.0),boundaryTags['obstacle'])]
    #add segments for square obstacle
    segments+=[(5,6,boundaryTags['obstacle']),
               (6,7,boundaryTags['obstacle']),
               (7,8,boundaryTags['obstacle']),
               (8,5,boundaryTags['obstacle'])]
    #obstacle
    points_on_obstacle=15
    radius = height/4.0
    center = (2.5,1.5)
# #     theta=0.0
#     pb  = cross_section(center,radius,theta)
#     vertices.append(('obstacle_'+'0',(pb[0],pb[1],boundaryTags['obstacle'])))
#     for gb in range(1,n_points_on_obstacle):
# #         theta = float(gb)/float(points_on_obstacle)*2.0*math.pi
#         pb = cross_secion(center,radius,theta)
#         vertices.append(('obstacle_'+`gb`,(pb[0],pb[1],boundaryTags['obstacle'])))
#         nv=len(vertices)    
#         segments.append((nv-2,nv-1,boundaryTags['obstacle']))
#     segments.append((nv-1,nv-points_on_obstacle,boundaryTags['obstacle']))
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
    poly.write("#vertices \n")
    for vN,v in enumerate(vertices):
        print (vN+1,v[1][0],v[1][1],v[2],v[0])
        poly.write('%d %12.5e %12.5e %d #%s \n' % (vN+1,v[1][0],v[1][1],v[2],v[0]))
    poly.write('%d %d \n' % (len(vertices),1))
    poly.write("#segments \n")
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0],s[1],s[2]))
    poly.write("#holes \n")
#     poly.write('%d\n' % (0,))
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          center[0],
                                          center[1],
                                          0+1))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    print vertices
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          vertices[0][1][0]+length*1.0e-8,
                                          vertices[0][1][1]+length*1.0e-8,
                                          0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('sophtest',cross_section=circular_cross_section)
    os.system('triangle -Aq30den sophtest.poly')
    os.system('showme sophtest.1.ele')
