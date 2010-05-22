#! /usr/bin/env python


#domain = [[0.0,4.0],[0.0,2.0]] #m

#hole   = [[1.5,2.5],[1.25,1.5]]


def genPoly(domain = [[0.0,4.0],[0.0,2.0]],#dimensions of box
            hole   = [[1.5,2.5],[1.25,1.5]],#and hole
            polyfileBase = 'hole2d'):
    
    X=0; Y=1
    L = (domain[X][1]-domain[X][0],domain[Y][1]-domain[Y][0],1.0)
    vertices = {}
    #vertices, starting at bottom going around, then inside box
    
    vertices[1] = (domain[X][0],domain[Y][0],0)
    vertices[2] = (domain[X][1],domain[Y][0],0)
    vertices[3] = (domain[X][1],domain[Y][1],0)
    vertices[4]= (domain[X][0],domain[Y][1],0)
    vertices[5]= (hole[X][0],hole[Y][0],2)
    vertices[6]= (hole[X][1],hole[Y][0],2)
    vertices[7]= (hole[X][1],hole[Y][1],2)
    vertices[8]= (hole[X][0],hole[Y][1],2)
    #
    boundaryTags = {'bottom':1,
                    'right':2,
                    'inflow':3,
                    'left':4,
                    'interior':0,
                    'hole bottom':5,
                    'hole top':5,
                    'hole right':5,
                    'hole left':5}

    #segments starting at bottom
    #seg id, endpoint 1, endpoint 2, id
    segments = {}
    segments[1]=(1,  2,  boundaryTags['bottom'], "#bottom ")
    segments[2]=(2,  3,  boundaryTags['right'], "#right ")
    segments[3]=(3,  4,  boundaryTags['inflow'], "#inflow ")
    segments[4]=(4,  1,  boundaryTags['left'], "#left ")
    segments[5]=(5,  6,  boundaryTags['hole bottom'], "#hole bottom ")
    segments[6]=(6,  7,  boundaryTags['hole right'], "#hole right ")
    segments[7]=(7,  8,  boundaryTags['hole top'], "#hole top ")
    segments[8]=(8,  5,  boundaryTags['hole left'], "#hole left ")

    #
    holes = {}
    holes[1] = (0.5*(hole[X][0]+hole[X][1]),0.5*(hole[Y][0]+hole[Y][1]))

    #regions
    regions = {}
    regions[1] = (0.5*(domain[X][0]+domain[X][1]),0.1*(hole[Y][0]-domain[Y][0])+domain[Y][0],0,"#background")


    nVertices = len(vertices)
    #base for keys should be one
    polyfile = open(polyfileBase+'.poly','w')
    polyfile.write('#poly file for 2d block hole example\n')
    polyfile.write('%d %d %d %d \n' % (nVertices,2,1,0))
    polyfile.write('#vertices \n')

    for iv,v in vertices.iteritems():
        polyfile.write('%d %12.5e %12.5e %d \n' % (iv,v[0],v[1],v[2]))

    nSegments = len(segments)
    polyfile.write('#segments \n')
    polyfile.write('%d %d \n' % (nSegments,1))

    for i,s in segments.iteritems():
        polyfile.write('%d %d %d %d %s \n' % (i,s[0],s[1],s[2],s[3]))

    nHoles = len(holes)
    polyfile.write('#holes \n')
    polyfile.write('%d \n' % (nHoles))
    for i,h in holes.iteritems():
        polyfile.write('%d %12.5e %12.5e\n' % (i,h[0],h[1]))

    nRegions = len(regions)
    polyfile.write('#regions \n')
    polyfile.write('%d \n' % (nRegions))
    for i,r in regions.iteritems():
        polyfile.write('%d %12.5e %12.5e %d %s \n' % (i,r[0],r[1],r[2],r[3]))

    polyfile.close()

    return L,boundaryTags

if __name__ == '__main__':
    domain = [[0.0,4.0],[0.0,2.0]] #m

    hole   = [[1.5,2.5],[1.5,1.75]]


    L,bt = genPoly(domain=domain,
                   hole=hole)

    print "domain= %s hole=%s L=%s boundaryTags=%s " % (domain,
                                                        hole,
                                                        L,
                                                        bt)

