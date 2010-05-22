import math

#This is the channel in plan view, rotated clockwise 90 degrees 
#            up
#           ____
#left_up   |    | right_up
#left_mid   )  (  right_mid
#left_down   ||   right_down
#            _
#           down
#y ->
#x
#|
#V

#couple of example profiles for curve
def linear_profile(x):
    return x
def sqrt_profile(x):
    return math.sqrt(x)
def quadratic_profile(x):
    return x**2
def cubic_profile(x):
    return -2.0*x**3+3.0*x**2
def cos_profile(x):
    return (1.0 - math.cos(x*math.pi))/2.0

#bottom must be planar for now
def flat_bottom(x,y):
    return 0.0

def genPoly(fileprefix,
            depth=1.0,
            up_width=3.0,
            up_to_down_ratio=3.0,
            up_len=1.0,
            mid_len=1.0,
            down_len=1.0,
            curve_fun=quadratic_profile,
            z_fun=flat_bottom,
            n_points_on_curve=10):
    boundaries=['upstream','downstream',
                'up_left','mid_left','down_left',
                'up_right','mid_right','down_right','top','bottom']
    boundaryTags=dict([(key,i) for (i,key) in enumerate(boundaries)])
    assert(curve_fun((0.0)) == 0.0)
    assert(curve_fun((1.0)) == 1.0)
    down_width = up_width/up_to_down_ratio
    down_x_start = up_len+mid_len
    down_y_start = up_width/2.0 - down_width/2.0
    flume_len = up_len+mid_len+down_len
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = {'upstream_left':(0,(0.0,0.0,z_fun(0.0,0.0)),boundaryTags['up_left']),
                'upstream_left_top':(1,(0.0,0.0,depth),boundaryTags['up_left']),
                'mid_start_left':(2,(up_len,0.0,z_fun(up_len,0.0)),boundaryTags['up_left']),
                'mid_start_left_top':(3,(up_len,0.0,depth),boundaryTags['up_left'])}
    facets = [ (0,2,3,1,boundaryTags['up_left']) ]
    #left curve
    dx = 1.0/float(n_points_on_curve-1)
    nv=4
    assert(nv==len(vertices))
    for i in range(1,n_points_on_curve-1): #go downstream
        xr_offset = i*dx
        x_offset = xr_offset*mid_len
        x = up_len+x_offset
        y = curve_fun(xr_offset)*down_y_start
        vertices['mid_left_'+`i`] = (nv,(x,y,z_fun(x,y)),boundaryTags['mid_left'])
        nv+=1
        vertices['mid_left_top_'+`i`]=(nv,(x,y,depth),boundaryTags['mid_left'])
        nv+=1
        assert(nv==len(vertices))
        facets.append( (nv-4,nv-2,nv-1,nv-3,boundaryTags['mid_left']))
    vertices.update({'mid_stop_left':(nv,(down_x_start,down_y_start,z_fun(down_x_start,down_y_start)),boundaryTags['down_left']),
                     'mid_stop_left_top':(nv+1,(down_x_start,down_y_start,depth),boundaryTags['down_left']),
                     'downstream_left':(nv+2,(flume_len,down_y_start,z_fun(flume_len,down_y_start)),boundaryTags['down_left']),
                     'downstream_left_top':(nv+3,(flume_len,down_y_start,depth),boundaryTags['down_left']),
                     'downstream_right':(nv+4,(flume_len,down_y_start+down_width,z_fun(flume_len,down_y_start+down_width)),boundaryTags['down_right']),
                     'downstream_right_top':(nv+5,(flume_len,down_y_start+down_width,depth),boundaryTags['down_right']),
                     'mid_stop_right':(nv+6,(down_x_start,down_y_start+down_width,z_fun(down_x_start,down_y_start+down_width)),boundaryTags['down_right']),
                     'mid_stop_right_top':(nv+7,(down_x_start,down_y_start+down_width,depth),boundaryTags['down_right'])})
    nv +=8
    assert(nv==len(vertices))
    facets += [(nv-10,nv-8,nv-7,nv-9,boundaryTags['mid_left']),
               (nv-8,nv-6,nv-5,nv-7,boundaryTags['down_left']),
               (nv-6,nv-4,nv-3,nv-5,boundaryTags['downstream']),
               (nv-4,nv-2,nv-1,nv-3,boundaryTags['down_right'])]
                

    #right curve (going back to upstream)
    for i in range(n_points_on_curve-2,0,-1): #go upstream
        xr_offset = i*dx
        x_offset = xr_offset*mid_len
        x = up_len+x_offset
        y = up_width-curve_fun(xr_offset)*(up_width-down_y_start-down_width)
        vertices['mid_right_'+`i`]=(nv,(x,y,z_fun(x,y)),
                                    boundaryTags['mid_right'])
        nv+=1
        vertices['mid_right_top_'+`i`] = (nv,(x,y,depth),
                                          boundaryTags['mid_right'])
        nv+=1
        assert(nv==len(vertices))
        facets.append( (nv-4,nv-2,nv-1,nv-3,boundaryTags['mid_right']))
    #
    vertices.update({'mid_start_right':(nv,(up_len,up_width,z_fun(up_len,up_width)),boundaryTags['up_right']),
                     'mid_start_right_top':(nv+1,(up_len,up_width,depth),boundaryTags['up_right']),
                     'upstream_right':(nv+2,(0.0,up_width,z_fun(0.0,up_width)),boundaryTags['up_right']),
                     'upstream_right_top':(nv+3,(0.0,up_width,depth),boundaryTags['up_right'])})
    nv+=4
    assert(nv==len(vertices))
    facets += [(nv-6,nv-4,nv-3,nv-5,boundaryTags['mid_right']),
               (nv-4,nv-2,nv-1,nv-3,boundaryTags['up_right']),
               (nv-2,0,1,nv-1,boundaryTags['upstream'])]
    bottomFacet = []
    topFacet=[]
    for vk,v in vertices.iteritems():
        if '_top' in vk:
            topFacet +=(v[0],)
        else:
            bottomFacet +=(v[0],)
    topFacet.sort()
    bottomFacet.sort()
    topFacet +=(boundaryTags['top'],)
    bottomFacet +=(boundaryTags['bottom'],)
    facets += [topFacet,bottomFacet]
    poly = open(fileprefix+'.poly','w')
    poly.write("#vertices \n")
    poly.write('%d %d %d %d \n' % (len(vertices),3,0,1))
    v_sorted = vertices.values()
    v_sorted.sort()
    for v in v_sorted:
        poly.write('%d %12.5e %12.5e %12.5e %d \n' % (v[0]+1,v[1][0],v[1][1],v[1][2],v[2]+1))
    poly.write("#facets \n")
    poly.write('%d %d \n' % (len(facets),1))
    for fN,f in enumerate(facets):
        poly.write('%d %d %d \n' % (1,0,f[-1]+1)) #1 polygon, 0 holes, tag
        nvf = len(f)-1
        facet_line = '%d' % (nvf,) #number  of vertices
        for i in range(nvf):
            facet_line += ' %d' % (f[i]+1,)
        facet_line+='\n'
        poly.write(facet_line)
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    #use a perturbation of the lower lefthand corner to get a pont inside
    poly.write('%d %12.5e %12.5e %12.5e %d\n' % (1,
                                                 vertices['upstream_left'][1][0]+1.0e-8*up_len,
                                                 vertices['upstream_left'][1][1]+1.0e-8*up_len,
                                                 vertices['upstream_left'][1][2]+1.0e-8*up_len,
                                                 0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('contraction_3d_linear',curve_fun=linear_profile)
    os.system('tetgen -Aq1.25en contraction_3d_linear.poly')
    os.system('tetview contraction_3d_linear.1.ele')
    
    genPoly('contraction_3d_linear_corner',mid_len=0.0,curve_fun=linear_profile)
    os.system('tetgen -Aq1.25en contraction_3d_linear_corner.poly')
    os.system('tetview contraction_3d_linear_corner.1.ele')
    
    genPoly('contraction_3d_sqrt',curve_fun=sqrt_profile)
    os.system('tetgen -Aq1.25en contraction_3d_sqrt.poly')
    os.system('tetview contraction_3d_sqrt.1.ele')

    genPoly('contraction_3d_quadratic',curve_fun=quadratic_profile)
    os.system('tetgen -Aq1.25en contraction_3d_quadratic.poly')
    os.system('tetview contraction_3d_quadratic.1.ele')
    
    genPoly('contraction_3d_cubic',curve_fun=cubic_profile)
    os.system(r'tetgen -Aq1.25en contraction_3d_cubic.poly')
    os.system(r'tetview contraction_3d_cubic.1.ele')

    genPoly('contraction_3d_cos',curve_fun=cos_profile)
    os.system(r'tetgen -Aq1.25en contraction_3d_cos.poly')
    os.system(r'tetview contraction_3d_cos.1.ele')

