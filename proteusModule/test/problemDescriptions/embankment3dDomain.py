#! /usr/bin/env python           
import math
from pyadh import Domain

def embankment3d(lengthBousDomain = 170.0,
                 ransDomainStop=1000.00,
                 ransDomainHeight = 10.0,
                 inflowLength = 10.0,
                 inflowPad  = 5.0,
                 outflowLength= 1.0,
                 outflowPad= 1.0,
                 width= 10.0):
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

    if inflowLength > inflowPad:
        print "Warning inflowLength= %s inflowPad= %s embankmentDomain setting inflowPad = %s " % (inflowLength,
                                                                                                   inflowPad,
                                                                                                   1.05*inflowLength)
        inflowPad = 1.05*inflowLength
    assert inflowLength < inflowPad
    #pop off points in Boussinesq domain
    bathymetryPoints = []
    for p in allBathymetryPoints:
        if p[0] >= lengthBousDomain and p[0] <= ransDomainStop:
            bathymetryPoints.append(p)

    backHeight = 4.54

    #have to assume points in correct order

    #pad for inflow outflow
    pin = bathymetryPoints[0]
    pout= bathymetryPoints[-1]

    #now get corners of domain
    pSW = [pin[0]-inflowPad,pin[1]]
    pSE = [pout[0]+outflowPad,pout[1]]
 
    ransDomainLength= pSE[0]-pSW[0]
    #domain
    L = (ransDomainLength,width,ransDomainHeight)
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

    #cek double check that all points are 2D, pSW -> [pSW[0],width,pSW[1]] ...

    frontVertices = [[p[0],0.0,p[1]] for p in bathymetryPoints]

    #start with NW corner and work way around
    #left 
    frontVertices.insert(0,[pNW[0],0.0,pNW[1]])
    frontVertices.insert(1,[pSW[0],0.0,pSW[1]])
    #add midpoint to make sure some points are inflow labelled
    #inflow is on bottom
    frontVertices.insert(2,[pSW[0]+0.5*inflowLength,0.0,pSW[1]])
    frontVertices.insert(3,[pSW[0]+inflowLength,0.0,pSW[1]])
    #
    frontVertices.append([bathymetryPoints[-1][0]+outflowPad-0.5*outflowLength,0.0,bathymetryPoints[-1][1]])
    #right
    frontVertices.append([pSE[0],0.0,pSE[1]])
    frontVertices.append([pNE[0],0.0,pNE[1]])
    backVertices = [[v[0],width,v[2]] for v in frontVertices]
    vertices = frontVertices+backVertices
    nfrontVertices = len(frontVertices)

    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    
    #cek build FACETS instead of segments. Each facet is a cyclic list of vertex, get rid of segment Flags and add one facet flag per facet
    facets = []
    facetFlags=[]
    front_facets=[[vN for vN in range(nfrontVertices)]]
    facets.append(front_facets)
    facetFlags.append(boundaryTags['front'])
    back_facets=[[vN for vN in range(nfrontVertices,2*nfrontVertices)]]
    facets.append(back_facets)
    print vertices
    facetFlags.append(boundaryTags['back'])
    top_facets=[[vN for vN in [0,nfrontVertices-1,nfrontVertices*2-1,nfrontVertices]]]#[[vN for vN in range(0,nfrontVertices+1)]]
    facets.append(top_facets)
    facetFlags.append(boundaryTags['top'])
    right_facets=[[vN for vN in [nfrontVertices-2,nfrontVertices-1,2*nfrontVertices-1,2*nfrontVertices-2]]]#[[vN for vN in range(1,1+(nfrontVertices+1))]]
    facets.append(right_facets)
    facetFlags.append(boundaryTags['right'])
    left_facets=[[vN for vN in [0,nfrontVertices,nfrontVertices+1,1]]]
    facets.append(left_facets)
    facetFlags.append(boundaryTags['left'])
    bottom_facets=[[vN,vN+1,vN+1+nfrontVertices,vN+nfrontVertices] for vN in range(1,nfrontVertices-1)]
    facets.append(bottom_facets)
    facetFlags.append(boundaryTags['bottom'])
    
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 facets=facets,
                                                 facetFlags=facetFlags)
    domain.backHeight = backHeight
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    domain =  embankment3d(lengthBousDomain=100.,ransDomainStop=230.0,ransDomainHeight=15.0,inflowLength=10.0)
    domain.writeAsymptote("embankment3d")
    domain.writePoly("embankment3d")
    os.system("asy -V embankment3d")
