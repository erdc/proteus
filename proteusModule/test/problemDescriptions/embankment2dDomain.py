#! /usr/bin/env python
import math
from pyadh import Domain

def embankment2d(lengthBousDomain = 170.0,
                 ransDomainStop=1000.00,
                 ransDomainHeight = 10.0,
                 inflowLength = 10.0,
                 inflowPad  = 5.0,
                 outflowLength= 1.0,
                 outflowPad= 1.0):
    """
    try to generate Jeff Melby's domain 
    points from sajlevee.bat file

    lengthBousDomain -- point in bathymetry data afterwhich RANS domain is considered
            
    ransDomainStop   -- point in bathymetry data afterwhich RANS domain stops

    inflowLength     -- size of region for inflow boundary
    inflowPad        -- how much to add to domain to accomodate inflow (make sure on flat region)
 
    outflowLength    -- size of region for outflow boundary
    outflowPad       -- how much to add to domain to accomodate outflow (make sure on flat region)
 
    
    """
    #lengthBousDomain = 170.0# where Boussinesq domain stops
    allBathymetryPoints = [[-100,        -8.53],
                           [101.71,	-8.53],
                           [120.00,	-7.92],
                           [132.19,	-7.52],
                           [184.01,	-5.79],
                           [190.41,	-3.66],
                           [196.81,	-1.52],
                           [204.43,	-1.52],
                           [204.43,	-1.22],
                           [205.34,	-1.22],
                           [205.34,	-0.91],
                           [206.26,	-0.91],
                           [206.26,	-0.61],
                           [207.17,	-0.61],
                           [207.17,	-0.30],
                           [208.09,	-0.30],
                           [208.09,	0.00],
                           [209.00,	0.00],
                           [209.00,	0.30],
                           [209.92,	0.30],
                           [209.92,	0.61],
                           [210.83,	0.61],
                           [210.83,	0.91],
                           [215.10,	0.91],
                           [215.10,	-1.00],
                           [250.00,	-1.25]]

    if inflowLength >= inflowPad:
        print "Warning inflowLength= %s inflowPad= %s embankmentDomain setting inflowPad = %s " % (inflowLength,
                                                                                                   inflowPad,
                                                                                                   1.05*inflowLength)
        inflowPad = 1.05*inflowLength
    assert inflowLength < inflowPad, "inflowLength= %s inflowPad=%s " % (inflowLength,inflowPad)
    #pop off points in Boussinesq domain
    bathymetryPoints = []
    for p in allBathymetryPoints:
        if p[0] >= lengthBousDomain and p[0] <= ransDomainStop:
            bathymetryPoints.append(p)

    backHeight = -1.25#4.54
    crestHeight=  0.91
    baseHeight = -1.52
    #have to assume points in correct order

    #pad for inflow outflow
    pin = bathymetryPoints[0]
    pout= bathymetryPoints[-1]

    #now get corners of domain
    pSW = [pin[0]-inflowPad,pin[1]]
    pSE = [pout[0]+outflowPad,pout[1]]
 
    ransDomainLength= pSE[0]-pSW[0]
    #domain
    L = (ransDomainLength,ransDomainHeight,1.0)
    #assume height of RANS domain gives enough vertical padding
    minY = min(pSW[1],pSE[1])
    pNW = [pSW[0],minY+ransDomainHeight]
    pNE = [pSE[0],minY+ransDomainHeight]

    
    #check vertical coordinates
    tmp = sorted(bathymetryPoints,cmp=lambda x,y: int(x[1]-y[1]))    
    assert minY <= tmp[0][1], "found point below proposed block floor minY=%s tmp[0]= " % (minY,tmp[0])
    assert minY+ransDomainHeight > tmp[-1][1], "found point above proposed block ceiling maxnY=%s tmp[-1]= " % (minY+ransDomainHeight,
                                                                                                                tmp[-1])
    #make our domain start at origin
    xshift=0.0-pSW[0]
    yshift=0.0-minY
    #backHeight += minY

    vertices = [p for p in bathymetryPoints]

    #start with NW corner and work way around
    #left 
    vertices.insert(0,pNW)
    vertices.insert(1,pSW)
    #add midpoint to make sure some points are inflow labelled
    #inflow is on bottom
    vertices.insert(2,[pSW[0]+0.5*inflowLength,pSW[1]])
    vertices.insert(3,[pSW[0]+inflowLength,pSW[1]])
    #
    vertices.append([bathymetryPoints[-1][0]+outflowPad-0.5*outflowLength,bathymetryPoints[-1][1]])
    #right
    vertices.append(pSE)
    vertices.append(pNE)
    nvertices = len(vertices)

    segmentLabels = {'left': 1,
                     'bottom' : 2,
                     'right'  : 3,
                     'top'    : 4,
                     'inflow' : 5,
                     'outflow': 6}

    segments = []
    segmentFlags=[]
    segments.append([0,1])
    segmentFlags.append(segmentLabels['left'])
    segments.append([1,2])
    segmentFlags.append(segmentLabels['inflow'])
    segments.append([2,3])
    segmentFlags.append(segmentLabels['inflow'])
    for i in range(3,nvertices-3):
        segments.append([i,i+1])
        segmentFlags.append(segmentLabels['bottom'])
    segments.append([nvertices-4,nvertices-3])
    segmentFlags.append(segmentLabels['outflow'])
    segments.append([nvertices-3,nvertices-2])
    segmentFlags.append(segmentLabels['outflow'])
    segments.append([nvertices-2,nvertices-1])
    segmentFlags.append(segmentLabels['right'])
    segments.append([nvertices-1,0])
    segmentFlags.append(segmentLabels['top'])
    print vertices,"s",segments,"sF",segmentFlags
    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags)
    domain.backHeight = backHeight
    domain.crestHeight= crestHeight
    domain.baseHeight = baseHeight
    #go ahead and add a boundary tags member 
    domain.boundaryTags = segmentLabels
    return domain

if __name__=='__main__':
    import os
    domain =  embankment2d(lengthBousDomain=170.,ransDomainStop=230.0,ransDomainHeight=10.,inflowLength=2.0)
    domain.writeAsymptote("embankment2d")
    domain.writePoly("embankment2d")
    os.system("asy -V embankment2d")
    print "embankment baseHeight= %s crestHeight= %s backHeight= %s " % (domain.baseHeight,domain.crestHeight,domain.backHeight)
