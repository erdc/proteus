#! /usr/bin/env python
import math


def genPoly(polyfileBase = "cobras_saj_embankment",
            lengthBousDomain = 170.0,ransDomainStop=1000.00,
            ransDomainHeight = 10.0,
            inflowLength = 10.0,
            inflowPad  = 15.0,
            outflowLength= 1.0,
            outflowPad= 1.0):
    """
    try to generate Jeff Melby's domain
    points from sajlevee.bat file

    lengthBousDomain -- point in bathymetry data afterwhich RANS domain is considered

    ransDomainStop   -- point in bathymetry data afterwhich RANS domain stops

    inflowLength     -- size of region for inflow boundary
    inflowPad        -- how much to add to domain to accomodate inflow (make sure on flat region)

    outflowLength    -- size of region for inflow boundary
    outflowPad       -- how much to add to domain to accomodate outflow (make sure on flat region)


    """
    #lengthBousDomain = 170.0# where Boussinesq domain stops
    allBathymetryPoints = [(-100,        -8.53),
                           (101.71,     -8.53),
                           (120.00,     -7.92),
                           (132.19,     -7.52),
                           (184.01,     -5.79),
                           (190.41,     -3.66),
                           (196.81,     -1.52),
                           (204.43,     -1.52),
                           (204.43,     -1.22),
                           (205.34,     -1.22),
                           (205.34,     -0.91),
                           (206.26,     -0.91),
                           (206.26,     -0.61),
                           (207.17,     -0.61),
                           (207.17,     -0.30),
                           (208.09,     -0.30),
                           (208.09,     0.00),
                           (209.00,     0.00),
                           (209.00,     0.30),
                           (209.92,     0.30),
                           (209.92,     0.61),
                           (210.83,     0.61),
                           (210.83,     0.91),
                           (215.10,     0.91),
                           (215.10,     -1.00),
                           (250.00,     -1.25)]


    #pop off points in Boussinesq domain
    bathymetryPoints = []
    for p in allBathymetryPoints:
        if p[0] >= lengthBousDomain and p[0] <= ransDomainStop:
            bathymetryPoints.append(p)

    backHeight = 4.54
    #how much to pad domain for inflow and outflow
    #inflowPad = 5.0 #m
    #outflowPad= 1.0
    #inflowLength = 1.0
    #outflowLength= 1.0

    #have to assume points in correct order

    #pad for inflow outflow
    pin = bathymetryPoints[0]
    pout= bathymetryPoints[-1]

    #now get corners of domain
    pSW = (pin[0]-inflowPad,pin[1])
    pSE = (pout[0]+outflowPad,pout[1])

    ransDomainLength= pSE[0]-pSW[0]
    #domain
    L = (ransDomainLength,ransDomainHeight,1.0)
    #assume height of RANS domain gives enough vertical padding
    minY = min(pSW[1],pSE[1])
    pNW = (pSW[0],minY+ransDomainHeight)
    pNE = (pSE[0],minY+ransDomainHeight)


    #check vertical coordinates
    tmp = sorted(bathymetryPoints,cmp=lambda x,y: int(x[1]-y[1]))
    assert minY <= tmp[0][1], "found point below proposed block floor minY=%s tmp[0]= " % (minY,tmp[0])
    assert minY+ransDomainHeight > tmp[-1][1], "found point above proposed block ceiling maxnY=%s tmp[-1]= " % (minY+ransDomainHeight,
                                                                                                                tmp[-1])

    #make our domain start at origin
    xshift=0.0-pSW[0]
    yshift=0.0-minY

    vertices = [p for p in bathymetryPoints]

    #start with NW corner and work way around
    #left
    vertices.insert(0,pNW)
    vertices.insert(1,pSW)
    #add midpoint to make sure some points are inflow labelled
    #inflow is on bottom
    vertices.insert(2,(pSW[0]+0.5*inflowLength,pSW[1]))
    vertices.insert(3,(pSW[0]+inflowLength,pSW[1]))

    #
    #vertices.append((bathymetryPoints[-1][0]+outflowPad-outflowLength,bathymetryPoints[-1][1]))
    vertices.append((bathymetryPoints[-1][0]+outflowPad-0.5*outflowLength,bathymetryPoints[-1][1]))
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
    segments.append((0,1,segmentLabels['left']))
    segments.append((1,2,segmentLabels['inflow']))
    segments.append((2,3,segmentLabels['inflow']))
    for i in range(3,nvertices-3):
        segments.append((i,i+1,segmentLabels['bottom']))
    segments.append((nvertices-4,nvertices-3,segmentLabels['outflow']))
    segments.append((nvertices-3,nvertices-2,segmentLabels['outflow']))
    segments.append((nvertices-2,nvertices-1,segmentLabels['right']))
    segments.append((nvertices-1,0,segmentLabels['top']))


    poly = open(polyfileBase+'.poly','w')
    poly.write('%d %d %d %d \n' % (nvertices,2,0,0))
    #write vertices
    poly.write("#vertices \n")
    for i,p in enumerate(vertices):
        poly.write('%d %12.5e %12.5e \n' % (i+1,xshift+p[0],yshift+p[1]))
    #write segments
    nSegments = len(segments)
    poly.write('%d %d \n' % (nSegments,1))
    poly.write("#segments \n")
    for sN,s in enumerate(segments):
        poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
    #if mesh just the outside of the structure insert holes here
    nholes = 0
    poly.write('%d \n' % (nholes,))
    poly.write("#holes \n")
    nsolid_regions = 0
    poly.write('%d \n' % (nholes,))
    poly.write("#solid regions \n")
    poly.close()

    return L,segmentLabels,backHeight
