import math
from pyadh import Domain

def pome2D( nx,
            ny,
            points_on_grain=50,
            points_on_boundary=50):
    """
    generate square lattice circle packing on unit square
    radius is determined by number of spheres
    returns boundary flags
    """
    #n_domain_vertices = 4
    #domain_vertices =[(0.0,0.0),(0.0,1.0),(1.0,1.0),(1.0,0.0)]
    radius = 1.0/(4.0*float(nx))

    points_on_arc= points_on_grain/2
    points_on_arc2=points_on_arc/4
    points_not = points_on_boundary/(nx)/2
    DA=3.0*radius/float(points_not)
    DB=2.0*radius/float(points_not)
    DC=1.0*radius/float(points_not/2)
    
    n_domain_vertices = 4 + 4*points_on_boundary
    DX=1.0/float(points_on_boundary+1)
    DY=1.0/float(points_on_boundary+1)
    domain_vertices=[]




    ####Top boundary
    for i in range(points_on_arc2+1):
        domain_vertices.append((radius*math.cos(math.pi*(1.0/4.0+float(i)/(4.0*float(points_on_arc2)))), radius*math.sin(math.pi*(1.0/4.0+float(i)/(4.0*float(points_on_arc2))))))
    for i in range(nx-1):
        for j in range(points_not):
            domain_vertices.append((0.0, domain_vertices[-1][1]+ DB))
        center = (0.0,4.0*radius+4.0*radius*float(i))
        for k in range(1,points_on_arc+1):    
            domain_vertices.append((center[0]+radius*math.cos(math.pi*(-1.0/2.0+1.0*float(k)/float(points_on_arc))), center[1]+radius*math.sin(math.pi*(-1.0/2.0+1.0*float(k)/float(points_on_arc)))))

    for j in range(points_not-1):
            domain_vertices.append((0.0, domain_vertices[-1][1]+ DB))
        #center = (0.0,4.0*radius+4.0*radius*float(i))
        #for k in range(1,points_on_arc+1):
    for i in range(points_on_arc2):
        domain_vertices.append((radius*math.cos(math.pi*(-1.0/2.0+float(i)/(4.0*float(points_on_arc2)))),1.0+ radius*math.sin(math.pi*(-1.0/2.0+float(i)/(4.0*float(points_on_arc2))))))

    L=len(domain_vertices)
    
    shift=(0.0,1.0)
    ang=math.pi/2.0
    for i in range(L):
        domain_vertices.append((shift[0]+domain_vertices[i][0]*math.cos(ang)+domain_vertices[i][1]*math.sin(ang),shift[1] -domain_vertices[i][0]*math.sin(ang)+domain_vertices[i][1]*math.cos(ang)))

    ang=math.pi
    shift=(1.0,1.0)

    for i in range(L):
       ##  domain_vertices.append((domain_vertices[i][0]*math.cos(ang)-domain_vertices[i][1]*math.sin(ang), domain_vertices[i][0]*math.sin(ang)+domain_vertices[i][1]*math.cos(ang))+shift)
        
        domain_vertices.append((shift[0]+domain_vertices[i][0]*math.cos(ang)+domain_vertices[i][1]*math.sin(ang),shift[1] -domain_vertices[i][0]*math.sin(ang)+domain_vertices[i][1]*math.cos(ang)))

    ang=3.0/2.0*math.pi
    shift=(1.0, 0.0)

    for i in range(L):
        ## domain_vertices.append((domain_vertices[i][0]*math.cos(ang)-domain_vertices[i][1]*math.sin(ang), domain_vertices[i][0]*math.sin(ang)+domain_vertices[i][1]*math.cos(ang))+shift)
        
        domain_vertices.append((shift[0]+domain_vertices[i][0]*math.cos(ang)+domain_vertices[i][1]*math.sin(ang),shift[1] -domain_vertices[i][0]*math.sin(ang)+domain_vertices[i][1]*math.cos(ang)))












    
 ##    ## for k in range(points_on_boundary+2):
## ##         domain_vertices.append((0.0,k*DY))
        
## ##     ## for k in range(1,points_on_boundary+2):
## ## ##         domain_vertices.append((k*DX,1.0))

## ## ##Top Boundary

## ##     point1=len(domain_vertices)
## ##     for i in range(1,points_not+1):
## ##         domain_vertices.append(( float(i)*DA,1.0))
## ##     for j in range(nx-1):
## ##         center=(4.0*radius+4.0*radius*float(j),1.0)
## ##         for k in range(1,points_on_arc+1):
## ##             domain_vertices.append((center[0]+radius*math.cos(math.pi*(1.0+float(k)/float(points_on_arc))), center[1]+radius*math.sin(math.pi*(1.0+float(k)/float(points_on_arc)))))
## ##         domain_vertices[-1]=(domain_vertices[-1][0],1.0)
## ##         for k in range(1,points_not+1):
## ##             domain_vertices.append((domain_vertices[-1][0]+DB,1.0))
## ##     for k in range(1,points_not/2 +1):
## ##         domain_vertices.append((domain_vertices[-1][0]+DC,1.0))

## ##     point2=len(domain_vertices)
## ##     domain_vertices[-1]=(1.0, 1.0)
    

## ##     # Right Boundary
## ##     for k in range(points_on_boundary,-1,-1):
## ##         domain_vertices.append((1.0,k*DY))

## ## #Bottom Boundary
## ##     for k in range(point2-point1-1):
## ##         domain_vertices.append((domain_vertices[point2-k-2][0], abs(-1.0*(domain_vertices[point2-k-2][1]-1.0))))
## ##    # domain_vertices.append((domain_vertices[point1][0], abs(-1.0*(domain_vertices[point1][1]-1.0))))

    
    
##     ## for i in range(1, points_not+1):
## ##         domain_vertices.append((domain_vertices[-1][0]-DA, 0.0))
## ##     for j in range(nx-1):
## ##         center=(1.0-4.0*radius-4.0*radius*float(j),0.0)
## ##         for k in range(1,points_on_arc+1):
## ##             domain_vertices.append((center[0]+radius*math.cos(math.pi*(1.0*float(k)/float(points_on_arc))), center[1]+radius*math.sin(math.pi*(1.0*float(k)/float(points_on_arc)))))
## ##         domain_vertices[-1]=(domain_vertices[-1][0], 0.0)
## ##         for k in range(1, points_not+1):
## ##             domain_vertices.append((domain_vertices[-1][0]-DB,0.0, 0.0))
## ##     for k in range(1,points_not/2 ):
## ##         domain_vertices.append((domain_vertices[-1][0]-DC,0.0)
#)
    ## for k in range(points_on_boundary,0,-1):
##         domain_vertices.append((k*DX,0.0))
    ## dx = 1.0/float(nx)
##     dy = 1.0/float(ny)
##     radius = 0.25*min(dx,dy)
    #radius = 1.0/(4.0*float(nx))
    grain_centers  = []
   
##     for i in range(ny):
##         for j in range(nx):
##             grain_centers.append((i*dy+0.5*dy,j*dx+0.5*dx))
    for i in range(nx):
        for j in range(nx):
                grain_centers.append((2.0*radius+float(i)*4.0*radius,2.0*radius+float(j)*4.0*radius))
    for i in range(nx-1):
        for j in range(nx-1):
            grain_centers.append((4.0*radius+float(i)*4.0*radius,4.0*radius+float(j)*4.0*radius)) 
            
            #grain_centers.append((2.0*radius+i*2.0*radius, 2.0*radius++j*2.0*radius))
        
    
    nvertices = len(domain_vertices) + len(grain_centers)*points_on_grain
    boundaries = ['left', 'right', 'front', 'back', 'obstacle']
    #based on vertex order and way initial segment list goes round
    boundaryFlags={'left':1,'back':2,'right':3,'front':4, 'obstacle':5}
    vertices=[]
    vertexFlags=[]
    segments=[]
    segmentFlags=[]
    #write vertices

    odds=[]
    for i in range(nx*2):
        odds.append(float(2*i+1)*radius)
    
    for v,p in enumerate(domain_vertices):#numbering is base 1 for triangle
        iscorner=False
        
        if (abs(p[0]) <1.0e-13):
            vertices.append([0.0,p[1]])
            for k in odds:
                if abs(k-p[1])<1.0e-5:
                    iscorner=True
            if iscorner==True:
                vertexFlags.append(boundaryFlags['obstacle'])
            else:
                vertexFlags.append(boundaryFlags['left'])
        elif (abs(p[0]- 1.0)<1.0e-13):
            vertices.append([1.0,p[1]])
            for k in odds:
                if abs(k-p[1])<1.0e-5:
                    iscorner=True
            if iscorner==True:
                vertexFlags.append(boundaryFlags['obstacle'])
            else:                  
                vertexFlags.append(boundaryFlags['right'])
        elif (abs(p[1]) <1.0e-13):
            vertices.append([p[0],0.0])
            for k in odds:
                if abs(k-p[0])<1.0e-5:
                    iscorner=True
            if iscorner==True:
                vertexFlags.append(boundaryFlags['obstacle'])
            else:
                vertexFlags.append(boundaryFlags['front'])
        elif (abs(p[1] - 1.0)<1.0e-13):
            #vertexFlags.append(boundaryFlags['back'])
            vertices.append([p[0],1.0])
            for k in odds:
                if abs(k-p[0])<1.0e-5:
                    iscorner=True
            if iscorner==True:
                vertexFlags.append(boundaryFlags['obstacle'])
            else:
                vertexFlags.append(boundaryFlags['back'])
        else:
            vertexFlags.append(boundaryFlags['obstacle'])
            vertices.append([p[0],p[1]])
    #write segments
    #left is X_minus, right is X_plus, front is Y_minus, back is Y_plus
    segments=[]
    segmentFlags=[]
    for sN in range(len(domain_vertices)-1):
        segments.append([sN,sN+1])
        if (abs(domain_vertices[sN][0]) < 1.0e-13 and abs(domain_vertices[sN+1][0]) <  1.0e-13):
            segmentFlags.append(boundaryFlags['left'])
        elif (abs(domain_vertices[sN][0] - 1.0) < 1.0e-13 and abs(domain_vertices[sN+1][0]-  1.0)< 1.0e-13):
            segmentFlags.append(boundaryFlags['right'])
        elif (abs(domain_vertices[sN][1])< 1.0e-13  and abs(domain_vertices[sN+1][1]) <  1.0e-13):
            segmentFlags.append(boundaryFlags['front'])
        elif (abs(domain_vertices[sN][1] - 1.0)<1.0e-13  and abs(domain_vertices[sN+1][1] -  1.0)<1.0e-13):
            segmentFlags.append(boundaryFlags['back'])
        elif (abs(domain_vertices[sN+1][0]) >1.0e-13 and abs(domain_vertices[sN+1][0] - 1.0)>1.0e-13):
            segmentFlags.append(boundaryFlags['obstacle'])
        elif (abs(domain_vertices[sN][0]) > 1.0e-13 and abs(domain_vertices[sN][0] - 1.0)>1.0e-13):
            segmentFlags.append(boundaryFlags['obstacle'])
        else:
            exit
    segments.append([len(domain_vertices)-1,0])
    segmentFlags.append(boundaryFlags['obstacle'])
   ##  if (domain_vertices[segments[-1][0]][0] == 0.0 and domain_vertices[segments[-1][1]][0] ==  0.0):
##         segmentFlags.append(boundaryFlags['left'])
##     if (domain_vertices[segments[-1][0]][0] == 1.0 and domain_vertices[segments[-1][1]][0] ==  1.0):
##         segmentFlags.append(boundaryFlags['right'])
##     if (domain_vertices[segments[-1][0]][1] == 0.0 and domain_vertices[segments[-1][1]][1] ==  0.0):
##         segmentFlags.append(boundaryFlags['front'])
##     if (domain_vertices[segments[-1][0]][1] == 1.0 and domain_vertices[segments[-1][1]][1] ==  1.0):
##         segmentFlags.append(boundaryFlags['back'])
    #end exterior boundary segments
    vStart = len(domain_vertices)
    sStart = len(segments)
    for g,c in enumerate(grain_centers):
        for gb in range(points_on_grain):
            vertices.append([c[0]+radius*math.sin(float(gb)/float(points_on_grain)*2.0*math.pi),c[1]+radius*math.cos(float(gb)/float(points_on_grain)*2.0*math.pi)])
            vertexFlags.append(boundaryFlags['obstacle'])
    for rb in range(len(grain_centers)):
        for gb in range(points_on_grain-1):
            segments.append([sStart+points_on_grain*rb+gb,sStart+points_on_grain*rb+gb+1])
            segmentFlags.append(boundaryFlags['obstacle'])
        segments.append([sStart+points_on_grain*rb+points_on_grain-1,sStart+points_on_grain*rb])
        segmentFlags.append(boundaryFlags['obstacle'])


    
    regions=[[vertices[0][0]+1.0e-8,
              vertices[0][1]+1.0e-8]]
    regionFlags=[1]

#construct domain object



    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  holes=grain_centers,
                                                  regions=regions,
                                                  regionFlags=regionFlags)
    #go ahead and add a boundary tags member 
    domain.boundaryFlags = boundaryFlags

    ## print len(vertices)
##     print len(segments)
##     print 'vertices = ' + str(vertices)
##     print 'vertex Flags' + str(vertexFlags)
##     print 'segments = '+ str(segments)
##     print 'segment Flags = ' +str(segmentFlags)

##     print len(vertices)
##     print len(vertexFlags)
##     print len(segments)
##     print len(segmentFlags)
              
    

    return domain

if __name__=='__main__':
    import os
    domain = pome2D(nx=7,ny=5, points_on_grain=35, points_on_boundary=50)
    domain.writeAsymptote("pome2D")
    domain.writePoly("pome2D")
    os.system("asy -V pome2D")
    
