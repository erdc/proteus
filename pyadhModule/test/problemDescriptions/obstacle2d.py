import math
from pyx import *
#This is the channel from the side
#              top
#             ____________
#            |            |
# upstream   | O          |  downstream
#            |            |
#             ____________
#             bottom
#^
#y x>


def circular_cross_section(center,radius,theta):
    return (radius*math.cos(theta)+center[0],
            radius*math.sin(theta)+center[1])

def genPoly(fileprefix,
            height=1.0,
            length=1.0,
            radius=0.25,
            center=(0.5,0.5),
            n_points_on_obstacle=10,
            cross_section=circular_cross_section):
    boundaries=['upstream',
                'downstream',
                'bottom',
                'top',
                'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [('upstream_bottom', (0.0,0.0),boundaryTags['bottom']),
                ('downstream_bottom',(length,0.0),boundaryTags['bottom']),
                ('downstream_top',(length,height),boundaryTags['top']),
                ('upstream_top',(0.0,height),boundaryTags['top'])]
    nv = len(vertices)
    segments = [(0,1,boundaryTags['bottom']),
                (1,2,boundaryTags['downstream']),
                (2,3,boundaryTags['top']),
                (3,0,boundaryTags['upstream'])]
    #cylinder
    theta=0.0
    pb  = cross_section(center,radius,theta)
    vertices.append(('obstacle_'+`0`,(pb[0],pb[1]),boundaryTags['obstacle']))
    for gb in range(1,n_points_on_obstacle):
        theta = float(gb)/float(n_points_on_obstacle)*2.0*math.pi
        pb  = cross_section(center,radius,theta)
        vertices.append(('obstacle_'+`gb`,(pb[0],pb[1]),boundaryTags['obstacle']))
        nv = len(vertices)
        segments.append((nv-2,nv-1,boundaryTags['obstacle']))
    segments.append((nv-1,nv-n_points_on_obstacle,boundaryTags['obstacle']))
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
    poly.write("#vertices \n")
    for vN,v in enumerate(vertices):
        poly.write('%d %12.5e %12.5e %d #%s \n' % (vN+1,v[1][0],v[1][1],v[2],v[0]))
    poly.write('%d %d \n' % (len(vertices),1))
    poly.write("#segments \n")
    g = graph.graphxy(width=8,
                      ratio=length/height,
                      x=graph.axis.linear(min=0,max=length),
                      y=graph.axis.linear(min=0,max=height))
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
        g.plot(graph.data.values(x=[vertices[s[0]][1][0],
                                    vertices[s[1]][1][0]],
                                 y=[vertices[s[0]][1][1],
                                    vertices[s[1]][1][1]]),
               styles=[graph.style.line()])
    g.writeEPSfile("obstacle2d")
    poly.write("#holes \n")
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e \n' % (1,
                                          center[0],
                                          center[1]))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e %d\n' % (1,
                                          vertices[0][1][0]+length*1.0e-8,
                                          vertices[0][1][1]+length*1.0e-8,
                                          0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
#     genPoly('cylinder_linear',cross_section=circular_cross_section)
#     os.system('triangle -Aq30den cylinder_linear.poly')
#     os.system('showme cylinder_linear.1.ele')
    genPoly('cylinder_test',height=2.5,length=5.0,radius=0.1,cross_section=circular_cross_section,center=(0.9,1.3))
    #os.system('triangle -Aq30den cylinder_test.poly')
    #os.system('showme cylinder_test.1.ele')
