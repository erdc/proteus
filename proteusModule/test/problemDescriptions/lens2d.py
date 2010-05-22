#! /usr/bin/env python


#domain = [[0.0,4.0],[0.0,2.0]] #m

#lens   = [[1.5,2.5],[1.25,1.5]]

#sourceXs = [1.89,1.9,2.1,2.11]

def genPoly(domain = [[0.0,4.0],[0.0,2.0]],#dimensions of box
            lens   = [[1.5,2.5],[1.25,1.5]],#and lens
            sourceXs = [1.89,1.9,2.1,2.11],#add extra locations to make sure get segment lables right
            polyfileBase = 'lens2d'):
    
    X=0; Y=1
    L = (domain[X][1]-domain[X][0],domain[Y][1]-domain[Y][0],1.0)
    vertices = {}
    #vertices, starting at bottom going around, then inside box
    
    vertices[1] = (domain[X][0],domain[Y][0],0)
    vertices[2] = (lens[X][0],domain[Y][0],0)
    vertices[3] = (lens[X][1],domain[Y][0],0)
    vertices[4] = (domain[X][1],domain[Y][0],0)
    vertices[5] = (domain[X][1],lens[Y][0],0)
    vertices[6] = (domain[X][1],lens[Y][1],0)
    vertices[7] = (domain[X][1],domain[Y][1],0)
    vertices[8] = (lens[X][1],domain[Y][1],0)
    vertices[9] = (sourceXs[-1],domain[Y][1],0)#pad source region to make sure gets correct label
    vertices[10]= (sourceXs[-2],domain[Y][1],1)
    vertices[11]= (sourceXs[-3],domain[Y][1],1)
    vertices[12]= (sourceXs[-4],domain[Y][1],0)
    vertices[13]= (lens[X][0],domain[Y][1],0)
    vertices[14]= (domain[X][0],domain[Y][1],0)
    vertices[15]= (domain[X][0],lens[Y][1],0)
    vertices[16]= (domain[X][0],lens[Y][0],0)
    vertices[17]= (lens[X][0],lens[Y][0],2)
    vertices[18]= (lens[X][1],lens[Y][0],2)
    vertices[19]= (lens[X][1],lens[Y][1],2)
    vertices[20]= (lens[X][0],lens[Y][1],2)
    #
    boundaryTags = {'bottom':1,
                    'right':2,
                    'top':3,
                    'inflow':4,
                    'left':5,
                    'interior':0,
                    'lens bottom':6,
                    'lens top':6,
                    'lens right':6,
                    'lens left':6}

    #segments starting at bottom
    #seg id, endpoint 1, endpoint 2, id
    segments = {}
    segments[1]=(1,  2,  boundaryTags['bottom'], "#bottom ")
    segments[2]=(2,  3,  boundaryTags['bottom'], "#bottom ")
    segments[3]=(3,  4,  boundaryTags['bottom'], "#bottom ")
    segments[4]=(4,  5,  boundaryTags['right'], "#right ")
    segments[5]=(5,  6,  boundaryTags['right'], "#right ")
    segments[6]=(6,  7,  boundaryTags['right'], "#right ")
    segments[7]=(7,  8,  boundaryTags['top'], "#top ")
    segments[8]=(8,  9,  boundaryTags['top'], "#top ")
    segments[9]=(9,  10, boundaryTags['inflow'], "#inflow pad ")
    segments[10]=(10, 11, boundaryTags['inflow'], "#inflow ")
    segments[11]=(11, 12, boundaryTags['inflow'], "#inflow pa ")
    segments[12]=(12, 13, boundaryTags['top'], "#top ")
    segments[13]=(13, 14, boundaryTags['top'], "#top ")
    segments[14]=(14, 15, boundaryTags['left'], "#left ")
    segments[15]=(15, 16, boundaryTags['left'], "#left ")
    segments[16]=(16, 1,  boundaryTags['left'], "#left, end outside ")
    segments[17]=(2, 17,  boundaryTags['interior'], "#outside inside ")
    segments[18]=(3,  18,  boundaryTags['interior'], "#outside inside ")
    segments[19]=(5,  18,  boundaryTags['interior'], "#outside inside ")
    segments[20]=(6,  19,  boundaryTags['interior'], "#outside inside ")
    segments[21]=(8,  19,  boundaryTags['interior'], "#outside inside ")
    segments[22]=(13, 20,  boundaryTags['interior'], "#outside inside ")
    segments[23]=(15, 20,  boundaryTags['interior'], "#outside inside ")
    segments[24]=(16, 17,  boundaryTags['interior'], "#outside inside ")
    segments[25]=(17, 18,  boundaryTags['lens bottom'], "#lens bottom")
    segments[26]=(18, 19,  boundaryTags['lens right'], "#lens right")
    segments[27]=(19, 20,  boundaryTags['lens top'], "#lens top")
    segments[28]=(20, 17,  boundaryTags['lens left'], "#lens left")

    #
    holes = {}

    #regions
    regions = {}
    regions[1] = (0.5*(domain[X][0]+domain[X][1]),0.1*(lens[Y][0]-domain[Y][0])+domain[Y][0],0,"#background")
    regions[2] = (0.5*(lens[X][0]+lens[X][1]),0.5*(lens[Y][0]+lens[Y][1]),1,"#lens")


    nVertices = len(vertices)
    #base for keys should be one
    polyfile = open(polyfileBase+'.poly','w')
    polyfile.write('#poly file for 2d block lense example\n')
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
        polyfile.write('%d %12.5e %12.5e %s \n' % (i,h[0],h[1],h[2]))

    nRegions = len(regions)
    polyfile.write('#regions \n')
    polyfile.write('%d \n' % (nRegions))
    for i,r in regions.iteritems():
        polyfile.write('%d %12.5e %12.5e %d %s \n' % (i,r[0],r[1],r[2],r[3]))

    polyfile.close()

    return L,boundaryTags

if __name__ == '__main__':
    domain = [[0.0,4.0],[0.0,2.0]] #m

    lens   = [[1.5,2.5],[1.25,1.5]]

    sourceXs = [1.89,1.9,2.1,2.11]

    L,bt = genPoly(domain=domain,
                   lens=lens,
                   sourceXs=sourceXs)

    print "domain= %s lens=%s sourceXs=%s L=%s boundaryTags=%s " % (domain,
                                                                    lens,
                                                                    sourceXs,
                                                                    L,
                                                                    bt)
    
