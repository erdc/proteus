import math

#This is the channel from the side, rotated clockwise 90 degrees 
#              up
#             ____
#            |    |
#bottom      | O  |  top
#            |    |
#             ____
#             down
#y ->
#x
#|
#V

def circular_cross_section(center,radius,theta):
    return (radius*math.sin(theta)+center[0],
            radius*math.cos(theta)+center[1])

def genPoly(fileprefix,
            height=1.0,
            length=1.0,
            width=1.0,
            radius=0.1,
            center=(0.5,0.5),
            n_points_on_obstacle=2*21-2,
            cross_section=circular_cross_section,
            thetaOffset=0.0):
    boundaries=['upstream',
                'downstream',
                'bottom',
                'top',
                'front',
                'back',
                'obstacle']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [('upstream_bottom', (0.0,0.0),boundaryTags['bottom']),
                ('downstream_bottom',(length,0.0),boundaryTags['bottom']),
                ('downstream_top',(length,height),boundaryTags['top']),
                ('upstream_top',(0.0,height),boundaryTags['top'])]
    nv = len(vertices)
    #cylinder
    theta=thetaOffset
    pb  = cross_section(center,radius,theta)
    vertices.append(('obstacle_'+`0`,(pb[0],pb[1]),boundaryTags['obstacle']))
    for gb in range(1,n_points_on_obstacle):
        theta = float(gb)/float(n_points_on_obstacle)*2.0*math.pi+thetaOffset
        pb  = cross_section(center,radius,theta)
        vertices.append(('obstacle_'+`gb`,(pb[0],pb[1]),boundaryTags['obstacle']))
    #now need to convert rep to 3D
    vertices3d={}
    front_cylinder=[]
    back_cylinder=[]
    for vN,v in enumerate(vertices):
        vertices3d[v[0]+'_front']=(vN,(v[1][0],0.0,v[1][1]),v[2])
        if 'obstacle' in v[0]:
            front_cylinder.append(vN)
    for vN,v in enumerate(vertices):
        vertices3d[v[0]+'_back']=(vN+len(vertices),(v[1][0],width,v[1][1]),v[2])
        if 'obstacle' in v[0]:
            back_cylinder.append(vN+len(vertices))
#     #try adding surface grid points
#     nnX=nnY=nnZ=5
#     dx = length/float(nnX-1)
#     dy = width/float(nnY-1)
#     dz = height/float(nnZ-1)
#     hardpoints=[]
#     for i in range(1,nnX-1):
#         for k in range(1,nnZ-1):
#             vertices3d['hardpoint'+`i*nnZ*nnY+k`]=(len(vertices3d)+1,(i*dx,0.0,k*dz),1)
#             vertices3d['hardpoint'+`i*nnZ*nnY+nnZ+k`]=(len(vertices3d)+1,(i*dx,width,k*dz),1)
#             hardpoints.append(vertices3d['hardpoint'+`i*nnZ*nnY+k`][0])
#             hardpoints.append(vertices3d['hardpoint'+`i*nnZ*nnY+nnZ+k`][0])
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices3d),3,0,1))
    poly.write("#vertices \n")
    a=vertices3d.values()
    a.sort()
    for v in a:
        poly.write('%d %12.5e %12.5e %12.5e %d #%s \n' % (1+v[0],v[1][0],v[1][1],v[1][2],v[2],v[-1]))
    poly.write("#facets \n")
#     poly.write('%d %d \n' % (6+n_points_on_obstacle+len(hardpoints),1))
    poly.write('%d %d \n' % (6+n_points_on_obstacle,1))
    #upstream
    poly.write('%d %d %d \n' % (1,0,boundaryTags['upstream'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_front'][0],
                                      1+vertices3d['upstream_bottom_back'][0],
                                      1+vertices3d['upstream_top_back'][0],
                                      1+vertices3d['upstream_top_front'][0]))
    #downstream
    poly.write('%d %d %d \n' % (1,0,boundaryTags['downstream'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['downstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_back'][0],
                                      1+vertices3d['downstream_top_back'][0],
                                      1+vertices3d['downstream_top_front'][0]))
    #top
    poly.write('%d %d %d \n' % (1,0,boundaryTags['top'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_top_front'][0],
                                      1+vertices3d['downstream_top_front'][0],
                                      1+vertices3d['downstream_top_back'][0],
                                      1+vertices3d['upstream_top_back'][0]))
    #bottom
    poly.write('%d %d %d \n' % (1,0,boundaryTags['bottom'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_back'][0],
                                      1+vertices3d['upstream_bottom_back'][0]))
    #front
    poly.write('%d %d %d \n' % (2,1,boundaryTags['front'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_front'][0],
                                      1+vertices3d['downstream_bottom_front'][0],
                                      1+vertices3d['downstream_top_front'][0],
                                      1+vertices3d['upstream_top_front'][0]))
    poly.write('%d' % n_points_on_obstacle)
    for vN in front_cylinder:
        poly.write(' %d' % (1+vN,))
    poly.write('\n')
    poly.write('%d %12.5e %12.5e %12.5e\n' % (1,center[0],0.0,center[1]))
    #back
    poly.write('%d %d %d \n' % (2,1,boundaryTags['back'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+vertices3d['upstream_bottom_back'][0],
                                      1+vertices3d['downstream_bottom_back'][0],
                                      1+vertices3d['downstream_top_back'][0],
                                      1+vertices3d['upstream_top_back'][0]))
    poly.write('%d' % n_points_on_obstacle)
    for vN in back_cylinder:
        poly.write(' %d' % (1+vN,))
    poly.write('\n')
    poly.write('%d %12.5e %12.5e %12.5e\n' % (1,center[0],width,center[1]))
    #cylinder
    for fN in range(n_points_on_obstacle-1):
        poly.write('%d %d %d \n' % (1,0,boundaryTags['obstacle'])) #1 polygon, 0 holes, tag
        poly.write('%d %d %d %d %d \n' % (4,
                                          1+front_cylinder[fN],
                                          1+back_cylinder[fN],
                                          1+back_cylinder[fN+1],
                                          1+front_cylinder[fN+1]))
    poly.write('%d %d %d \n' % (1,0,boundaryTags['obstacle'])) #1 polygon, 0 holes, tag
    poly.write('%d %d %d %d %d \n' % (4,
                                      1+front_cylinder[-1],
                                      1+back_cylinder[-1],
                                      1+back_cylinder[0],
                                      1+front_cylinder[0]))
    #hardpoints
#     for p in hardpoints:
#         poly.write('%d %d %d \n' % (1,0,1)) #1 polygon, 0 holes, tag
#         poly.write('%d %d \n' % (1,
#                                  p))
    #finished with facets
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    #use a perturbation of the lower bottomhand corner to get a pont inside
    poly.write('%d %12.5e %12.5e %12.5e %d\n' % (1,
                                                 1.0e-8*height,
                                                 1.0e-8*height,
                                                 1.0e-8*height,
                                                 0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('cylinder_3d_linear',n_points_on_obstacle=4,thetaOffset=math.pi/4.0)
    os.system('tetgen -YAq1.25en cylinder_3d_linear.poly')
    os.system('tetview cylinder_3d_linear.1.ele')
