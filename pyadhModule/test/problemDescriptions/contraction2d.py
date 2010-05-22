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

def genPoly(fileprefix,
            up_width=3.0,
            up_to_down_ratio=3.0,
            up_len=1.0,
            mid_len=1.0,
            down_len=1.0,
            curve_fun=quadratic_profile,
            n_points_on_curve=10):
    boundaries=['upstream','downstream',
                'up_left','mid_left','down_left',
                'up_right','mid_right','down_right']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    assert(curve_fun((0.0)) == 0.0)
    assert(curve_fun((1.0)) == 1.0)
    down_width = up_width/up_to_down_ratio
    down_x_start = up_len+mid_len
    down_y_start = up_width/2.0 - down_width/2.0
    flume_len = up_len+mid_len+down_len
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [('upstream_left', (0.0,0.0),boundaryTags['up_left']),
                ('mid_start_left',(up_len,0.0),boundaryTags['up_left'])]
    segments = [ (0,1,boundaryTags['up_left']) ]
    #left curve
    dx = 1.0/float(n_points_on_curve-1)
    for i in range(1,n_points_on_curve): #go downstream
        xr_offset = i*dx
        x_offset = xr_offset*mid_len
        vertices.append(('mid_left_'+`i`,(up_len+x_offset,
                                          curve_fun(xr_offset)*down_y_start),boundaryTags['mid_left']))
        nv = len(vertices)
        segments.append( (nv-2,nv-1,boundaryTags['mid_left']))
    vertices += [('mid_stop_left',(down_x_start,down_y_start),boundaryTags['down_left']),
                 ('downstream_left',(flume_len,down_y_start),boundaryTags['down_left']),
                 ('downstream_right',(flume_len,down_y_start+down_width),boundaryTags['down_right']),
                 ('mid_stop_right',(down_x_start,down_y_start+down_width),boundaryTags['down_right'])]
    nv = len(vertices)
    segments += [(nv-5,nv-4,boundaryTags['mid_left']),
                 (nv-4,nv-3,boundaryTags['down_left']),
                 (nv-3,nv-2,boundaryTags['downstream']),
                 (nv-2,nv-1,boundaryTags['down_right'])]
    #right curve (going back to upstream)
    for i in range(n_points_on_curve-1,0,-1): #go upstream
        xr_offset = i*dx
        x_offset = xr_offset*mid_len
        vertices.append(('mid_right_'+`i`,(up_len+x_offset,
                                           up_width-curve_fun(xr_offset)*(up_width-down_y_start-down_width)),
                         boundaryTags['mid_right']))
        nv = len(vertices)
        segments.append( (nv-2,nv-1,boundaryTags['mid_right']))
    #
    vertices += [('mid_start_right',(up_len,up_width),boundaryTags['up_right']),
                 ('upstream_right',(0.0,up_width),boundaryTags['up_right'])]
    nv = len(vertices)
    segments += [(nv-3,nv-2,boundaryTags['mid_right']),
                 (nv-2,nv-1,boundaryTags['up_right']),
                 (nv-1,0,boundaryTags['upstream'])]
    poly = open(fileprefix+'.poly','w')
    poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
    poly.write("#vertices \n")
    for vN,v in enumerate(vertices):
        poly.write('%d %12.5e %12.5e %d #%s \n' % (vN+1,v[1][0],v[1][1],v[2],v[0]))
    poly.write('%d %d \n' % (len(vertices),1))
    poly.write("#segments \n")
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
    poly.write("#segments \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e %d' % (1,
                                        vertices[0][1][0]+up_len*1.0e-8,
                                        vertices[0][1][1]+up_len*1.0e-8,
                                        0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('contraction_linear',curve_fun=linear_profile)
    os.system('triangle -Aq30den contraction_linear.poly')
    os.system('showme contraction_linear.1.ele')
    
    genPoly('contraction_linear_corner',mid_len=0.0,curve_fun=linear_profile)
    os.system('triangle -Aq30den contraction_linear_corner.poly')
    os.system('showme contraction_linear_corner.1.ele')
    
    genPoly('contraction_sqrt',curve_fun=sqrt_profile)
    os.system('triangle -Aq30den contraction_sqrt.poly')
    os.system('showme contraction_sqrt.1.ele')

    genPoly('contraction_quadratic',curve_fun=quadratic_profile)
    os.system('triangle -Aq30den contraction_quadratic.poly')
    os.system('showme contraction_quadratic.1.ele')
    
    genPoly('contraction_cubic',curve_fun=cubic_profile)
    os.system(r'triangle -Aq30den contraction_cubic.poly')
    os.system(r'showme contraction_cubic.1.ele')
    
    genPoly('contraction_cos',curve_fun=cos_profile)
    os.system(r'triangle -Aq30den contraction_cos.poly')
    os.system(r'showme contraction_cos.1.ele')

