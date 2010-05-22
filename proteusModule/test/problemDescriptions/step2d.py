import math

#This is the channel from the side, rotated clockwise 90 degrees 
#
#                  top
#        ___________________
#        |                  |
#upstream|                  |
#        _______            | downstream
#upstream bottom|           |
#           step|           |
#               |___________|
#          downstream bottom
#^
#y x>

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
            upstream_height=0.5,
            downstream_height=1.0,
            upstream_length=1.0,
            downstream_length=5.0,
            step_fun=quadratic_profile,
            n_points_on_step=10,
            step_length=1.0):
    boundaries=['upstream',
                'downstream',
                'upstream_bottom',
                'downstream_bottom',
                'top',
                'step_bottom']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    assert(step_fun((0.0)) == 0.0)
    assert(step_fun((1.0)) == 1.0)
    length=upstream_length+step_length+downstream_length
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = [('upstream_bottom', (0.0,downstream_height-upstream_height),boundaryTags['upstream_bottom']),
                ('step_bottom',(upstream_length,downstream_height-upstream_height),boundaryTags['upstream_bottom'])]
    segments = [ (0,1,boundaryTags['upstream_bottom']) ]
    #bottom curve
    dx = 1.0/float(n_points_on_step-1)
    for i in range(1,n_points_on_step-1): #go downstream
        xr_offset = i*dx
        x_offset = xr_offset*step_length
        vertices.append(('step_'+`i`,(upstream_length+x_offset,
                                      step_fun(1.0-xr_offset)*(downstream_height-upstream_height)),boundaryTags['step_bottom']))
        nv = len(vertices)
        segments.append( (nv-2,nv-1,boundaryTags['step_bottom']))
    vertices += [('step_'+`n_points_on_step`,(upstream_length+step_length,0.0),boundaryTags['downstream_bottom']),
                 ('downstream_bottom',(length,0.0),boundaryTags['downstream_bottom']),
                 ('downstream_top',(length,downstream_height),boundaryTags['top']),
                 ('upstream_top',(0.0,downstream_height),boundaryTags['top'])]
    nv = len(vertices)
    segments += [(nv-5,nv-4,boundaryTags['step_bottom']),
                 (nv-4,nv-3,boundaryTags['downstream_bottom']),
                 (nv-3,nv-2,boundaryTags['downstream']),
                 (nv-2,nv-1,boundaryTags['top']),
                 (nv-1,0,boundaryTags['upstream'])]
    poly = open(fileprefix+'.poly','w')
    poly.write("#vertices \n")
    poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
    for vN,v in enumerate(vertices):
        poly.write('%d %12.5e %12.5e %d #%s \n' % (vN+1,v[1][0],v[1][1],v[2],v[0]))
    poly.write("#segments \n")
    poly.write('%d %d \n' % (len(segments),1))
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
    poly.write("#holes \n")
    poly.write('%d\n' % (0,))
    poly.write("#regions \n")
    poly.write('%d\n' % (1,))
    poly.write('%d %12.5e %12.5e %d' % (1,
                                        vertices[0][1][0]+upstream_length*1.0e-8,
                                        vertices[0][1][1]+upstream_length*1.0e-8,
                                        0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    genPoly('step_linear',step_fun=linear_profile)
    os.system('triangle -Aq30den step_linear.poly')
    os.system('showme step_linear.1.ele')
    
    genPoly('step_linear_corner',
            step_fun=linear_profile,
            n_points_on_step=2,
            step_length=0.0)
    os.system('triangle -Aq30den step_linear_corner.poly')
    os.system('showme step_linear_corner.1.ele')
    
    genPoly('step_sqrt',step_fun=sqrt_profile)
    os.system('triangle -Aq30den step_sqrt.poly')
    os.system('showme step_sqrt.1.ele')

    genPoly('step_quadratic',step_fun=quadratic_profile)
    os.system('triangle -Aq30den step_quadratic.poly')
    os.system('showme step_quadratic.1.ele')
    
    genPoly('step_cubic',step_fun=cubic_profile)
    os.system(r'triangle -Aq30den step_cubic.poly')
    os.system(r'showme step_cubic.1.ele')
    
    genPoly('step_cos',step_fun=cos_profile)
    os.system(r'triangle -Aq30den step_cos.poly')
    os.system(r'showme step_cos.1.ele')

