import math
from pyadh import Domain

def ellipse2D( nx,
               ny,
               b_over_a,
               tilt_angle,
               points_on_grain,
               points_on_boundary):
    """
    generate square lattice ellipse packing on unit square
    axis lengths are determined by number of spheres
    proportion of axes are given by b_over_a
    returns boundary flags
    """
    #n_domain_vertices = 4
    #domain_vertices =[(0.0,0.0),(0.0,1.0),(1.0,1.0),(1.0,0.0)]
    n_domain_vertices = 4 + 4*points_on_boundary
    DX=1.0/float(points_on_boundary+1)
    DY=1.0/float(points_on_boundary+1)
    domain_vertices=[]
    for k in range(points_on_boundary+2):
        domain_vertices.append((0.0,k*DY))
    for k in range(1,points_on_boundary+2):
        domain_vertices.append((k*DX,1.0))
    for k in range(points_on_boundary,-1,-1):
        domain_vertices.append((1.0,k*DY))
    for k in range(points_on_boundary,0,-1):
        domain_vertices.append((k*DX,0.0))
    dx = 1.0/float(nx)
    dy = 1.0/float(ny)
#    radius = 0.25*min(dx,dy)
    
    xax = 0.25*min(dx,dy)
    yax = xax*b_over_a
    
    grain_centers  = []
   
    for i in range(ny):
        for j in range(nx):
            grain_centers.append((i*dy+0.5*dy,j*dx+0.5*dx))
    nvertices = len(domain_vertices) + len(grain_centers)*points_on_grain
    boundaries = ['left', 'right', 'front', 'back', 'obstacle']
    #based on vertex order and way initial segment list goes round
    boundaryFlags={'left':1,'back':2,'right':3,'front':4, 'obstacle':5}
    vertices=[]
    vertexFlags=[]
    segments=[]
    segmentFlags=[]
    #write vertices
    
    for v,p in enumerate(domain_vertices):#numbering is base 1 for triangle
        if (p[0] == 0.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['left'])
        elif (p[0] == 1.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['right'])
        elif (p[1] == 0.0):
            vertices.append([p[0],p[1]])
            vertexFlags.append(boundaryFlags['front'])
        elif (p[1] == 1.0):
            vertexFlags.append(boundaryFlags['back'])
            vertices.append([p[0],p[1]])
        else:
            exit
    #write segments
    #left is X_minus, right is X_plus, front is Y_minus, back is Y_plus
    segments=[]
    segmentFlags=[]
    for sN in range(len(domain_vertices)-1):
        segments.append([sN,sN+1])
        if (domain_vertices[sN][0] == 0.0 and domain_vertices[sN+1][0] ==  0.0):
            segmentFlags.append(boundaryFlags['left'])
        elif (domain_vertices[sN][0] == 1.0 and domain_vertices[sN+1][0] ==  1.0):
            segmentFlags.append(boundaryFlags['right'])
        elif (domain_vertices[sN][1] == 0.0 and domain_vertices[sN+1][1] ==  0.0):
            segmentFlags.append(boundaryFlags['front'])
        elif (domain_vertices[sN][1] == 1.0 and domain_vertices[sN+1][1] ==  1.0):
            segmentFlags.append(boundaryFlags['back'])
        else:
            exit
    segments.append([len(domain_vertices)-1,0])
    if (domain_vertices[segments[-1][0]][0] == 0.0 and domain_vertices[segments[-1][1]][0] ==  0.0):
        segmentFlags.append(boundaryFlags['left'])
    if (domain_vertices[segments[-1][0]][0] == 1.0 and domain_vertices[segments[-1][1]][0] ==  1.0):
        segmentFlags.append(boundaryFlags['right'])
    if (domain_vertices[segments[-1][0]][1] == 0.0 and domain_vertices[segments[-1][1]][1] ==  0.0):
        segmentFlags.append(boundaryFlags['front'])
    if (domain_vertices[segments[-1][0]][1] == 1.0 and domain_vertices[segments[-1][1]][1] ==  1.0):
        segmentFlags.append(boundaryFlags['back'])
    #end exterior boundary segments
    vStart = len(domain_vertices)
    sStart = len(segments)
    for g,c in enumerate(grain_centers):
        for gb in range(points_on_grain):
            a1= xax*math.sin(float(gb)/float(points_on_grain)*2.0*math.pi)
            a2= yax*math.cos(float(gb)/float(points_on_grain)*2.0*math.pi)
            vertices.append([c[0]+a1*math.cos(tilt_angle) - a2*math.sin(tilt_angle),
                             c[1]-a1*math.sin(tilt_angle) + a2*math.cos(tilt_angle)])                 
           # vertices.append([c[0]+xax*math.sin(float(gb)/float(points_on_grain)*2.0*math.pi),c[1]+yax*math.cos(float(gb)/float(points_on_grain)*2.0*math.pi)])
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

    return domain

if __name__=='__main__':
    import os
    domain = ellipse2D(nx=5,ny=5,b_over_a =1.0/0.80,tilt_angle= 0.0, points_on_grain=50, points_on_boundary=50)
    domain.writeAsymptote("ellipse2D")
    domain.writePoly("ellipse2D")
    os.system("asy -V ellipse2D")
