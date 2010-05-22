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

#couple of example profiles for step
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
            width=0.5,
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
                'step_bottom',
                'downstream_bottom',
                'top',
                'back',
                'front']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    assert(step_fun((0.0)) == 0.0)
    assert(step_fun((1.0)) == 1.0)
    length =upstream_length+step_length+downstream_length
    upstream_bottom = downstream_height-upstream_height
    dx = 1.0/float(n_points_on_step-1)
    #work around the domain from (0.0,0.0) going counterclockwise
    vertices = {'upstream_bottom':(0,(0.0,0.0,upstream_bottom),boundaryTags['upstream_bottom']),
                'upstream_bottom_back':(1,(0.0,width,upstream_bottom),boundaryTags['upstream_bottom']),
                'step_start_bottom':(2,(upstream_length,0,upstream_bottom),boundaryTags['upstream_bottom']),
                'step_start_bottom_back':(3,(upstream_length,width,upstream_bottom),boundaryTags['upstream_bottom'])}
    facets = [ (0,2,3,1,boundaryTags['upstream_bottom']) ]
    #bottom curve
    dx = 1.0/float(n_points_on_step-1)
    nv=4
    assert(nv==len(vertices))
    for i in range(1,n_points_on_step-1): #go downstream
        xr_offset = i*dx
        x_offset = xr_offset*step_length
        x = upstream_length+x_offset
        vertices['step_bottom'+`i`] = (nv,(x,0.0,step_fun(1.0-xr_offset)*(upstream_bottom)),boundaryTags['step_bottom'])
        nv+=1
        vertices['step_bottom_back_'+`i`]=(nv,(x,width,step_fun(1.0-xr_offset)*(upstream_bottom)),boundaryTags['step_bottom'])
        nv+=1
        assert(nv==len(vertices))
        facets.append( (nv-4,nv-2,nv-1,nv-3,boundaryTags['step_bottom']))
    vertices.update({'step_bottom_stop':(nv,(upstream_length+step_length,0.0,0.0),boundaryTags['downstream_bottom']),
                     'step_bottom_stop_back':(nv+1,(upstream_length+step_length,width,0.0),boundaryTags['downstream_bottom']),
                     'downstream_bottom':(nv+2,(length,0.0,0.0),boundaryTags['downstream_bottom']),
                     'downstream_bottom_back':(nv+3,(length,width,0.0),boundaryTags['downstream_bottom']),
                     'downstream_top':(nv+4,(length,0.0,downstream_height),boundaryTags['top']),
                     'downstream_top_back':(nv+5,(length,width,downstream_height),boundaryTags['top']),
                     'upstream_top':(nv+6,(0.0,0.0,downstream_height),boundaryTags['top']),
                     'upstream_top_back':(nv+7,(0.0,width,downstream_height),boundaryTags['top'])})
    nv +=8
    assert(nv==len(vertices))
    facets += [(nv-10,nv-8,nv-7,nv-9,boundaryTags['step_bottom']),
               (nv-8,nv-6,nv-5,nv-7,boundaryTags['downstream_bottom']),
               (nv-6,nv-4,nv-3,nv-5,boundaryTags['downstream']),
               (nv-4,nv-2,nv-1,nv-3,boundaryTags['top']),
               (nv-2,0,1,nv-1,boundaryTags['upstream'])]
    frontFacet = []
    backFacet=[]
    for vk,v in vertices.iteritems():
        if '_back' in vk:
            backFacet +=(v[0],)
        else:
            frontFacet +=(v[0],)
    backFacet.sort()
    frontFacet.sort()
    backFacet +=(boundaryTags['back'],)
    frontFacet +=(boundaryTags['front'],)
    facets += [backFacet,frontFacet]
    poly = open(fileprefix+'.poly','w')
    poly.write("#vertices \n")
    poly.write('%d %d %d %d \n' % (len(vertices),3,0,1))
    v_sorted = vertices.values()
    v_sorted.sort()
    for v in v_sorted:
        poly.write('%d %12.5e %12.5e %12.5e %d \n' % (v[0]+1,v[1][0],v[1][1],v[1][2],v[2]))
    poly.write("#facets \n")
    poly.write('%d %d \n' % (len(facets),1))
    for fN,f in enumerate(facets):
        poly.write('%d %d %d \n' % (1,0,f[-1])) #1 polygon, 0 holes, tag
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
    #use a perturbation of the lower bottomhand corner to get a pont inside
    poly.write('%d %12.5e %12.5e %12.5e %d\n' % (1,
                                                 vertices['upstream_bottom'][1][0]+1.0e-8*upstream_length,
                                                 vertices['upstream_bottom'][1][1]+1.0e-8*upstream_length,
                                                 vertices['upstream_bottom'][1][2]+1.0e-8*upstream_length,
                                                 0+1))
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    upstream_height=0.5
    downstream_height=1.0
    upstream_length = 1.0
    downstream_length = 5.0
    length = upstream_length+downstream_length
    polyfile = "step3d_test"
    triangleOptions = "VpAq1.25Dena%e" % ((0.15*upstream_height)**3 / 6.0)
    boundaryTags = genPoly(fileprefix=polyfile,
                           width=upstream_height/10.0,
                           upstream_height=upstream_height,
                           downstream_height=downstream_height,
                           upstream_length=upstream_length,
                           downstream_length=downstream_length,
                           step_fun=linear_profile,
                           n_points_on_step=2,
                           step_length=0.0)
    os.system('tetgen -%s step3d_test.poly' % (triangleOptions,))
    os.system('tetview step3d_test.1.ele')
    
#     genPoly('step_3d_linear_corner',step_length=0.0,step_fun=linear_profile,n_points_on_step=2)
#     os.system('tetgen -YAq1.25en step_3d_linear_corner.poly')
#     os.system('tetview step_3d_linear_corner.1.ele')
    
#     genPoly('step_3d_sqrt',step_fun=sqrt_profile)
#     os.system('tetgen -YAq1.25en step_3d_sqrt.poly')
#     os.system('tetview step_3d_sqrt.1.ele')

#     genPoly('step_3d_quadratic',step_fun=quadratic_profile)
#     os.system('tetgen -YAq1.25en step_3d_quadratic.poly')
#     os.system('tetview step_3d_quadratic.1.ele')
    
#     genPoly('step_3d_cubic',step_fun=cubic_profile)
#     os.system(r'tetgen -YAq1.25en step_3d_cubic.poly')
#     os.system(r'tetview step_3d_cubic.1.ele')

#     genPoly('step_3d_cos',step_fun=cos_profile)
#     os.system(r'tetgen -YAq1.25en step_3d_cos.poly')
#     os.system(r'tetview step_3d_cos.1.ele')

