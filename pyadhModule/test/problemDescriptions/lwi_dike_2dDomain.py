#! /usr/bin/env python
import math
from pyadh import Domain

"""
 generate waves and runup on LWI dike from OPTICREST report

 units should be meters by default


 z=0         beach slope            back slope  
             (1/6)            B      (-1/3)
 |--------|----------------|-------|---------|
 0        x_bs            x_be    x_ce       x_bse=L

x_bs = inflow length (1 [m] default), 
z(x) = 0, [0,x_bs]

x_be = x_bs + beachLength (4.8 [m])
z(x) = 0.0 + (x-x_bs)beachSlope, [x_bs,x_be]

x_ce = x_be + B [m]
z(x) = z(x_be), [x_be,x_ce]

x_bse= x_ce + backSlopeLength
z(x) = z(x_ce) + (x-x_ce)*backSlope [x_ce,x_bse]

x_L  = x_bse [m]
z(x) = z(x_bse)

total domain height = z(x_bse) + domainHeightPad
"""

def lwi_dike_2d(domainHeightPad=0.27,
                inflowLength=1.0,
                beachLength=4.8,
                beachSlope =0.167,
                B  = 0.3,
                backSlopeLength=0.2,
                backSlope=-0.333):

    """
    """
    #describe bathymetry as above
    x_0  = 0.0
    x_bs = x_0 + inflowLength; x_be = x_bs + beachLength
    x_ce = x_be + B
    x_bse= x_ce + backSlopeLength

    def bathymetry(x):
        if x[0] <= x_bs:  return 0.0
        if x[0] <= x_be:  return 0.0 + (x[0]-x_bs)*beachSlope
        if x[0] <= x_ce:  return bathymetry([x_be])
        if x[0] <= x_bse: return bathymetry([x_ce]) + (x[0]-x_ce)*backSlope
        return bathymetry([x_bse])
    def bathymetryGrad(x):
        if x[0] <= x_bs:  return (0.0,)
        if x[0] <= x_be:  return (beachSlope,) #beach slope
        if x[0] <= x_ce:  return (0.0,)
        if x[0] <= x_bse: return (backSlope,)
        return (0.0,)
    
    x_L  = x_bse
    z_L  = bathymetry([x_bse]) + domainHeightPad
    L = [x_L,z_L]

    xpoints = [x_0,x_bs,x_be,x_ce,x_bse]
    zpoints = [bathymetry([x]) for x in xpoints]
    bathymetryPoints = [[x,z] for x,z in zip(xpoints,zpoints)]

    backHeight = bathymetry([x_bse])

    #have to assume points in correct order

    #pad for inflow outflow
    pSW = bathymetryPoints[0]
    pSE= bathymetryPoints[-1]

    #now get top corners of domain
    minY = 0.0
    pNW = [pSW[0],minY+z_L]
    pNE = [pSE[0],minY+z_L]

    
    #check vertical coordinates
    tmp = sorted(bathymetryPoints,cmp=lambda x,y: int(x[1]-y[1]))    
    assert minY <= tmp[0][1], "found point below proposed block floor minY=%s tmp[0]= " % (minY,tmp[0])
    assert minY+ z_L > tmp[-1][1], "found point above proposed block ceiling maxnY=%s tmp[-1]= " % (minY+z_L,
                                                                                                    tmp[-1])
    vertices = [p for p in bathymetryPoints]

    #start with NW corner and work way around
    #left 
    vertices.insert(0,pNW)
    #not needed if no pad vertices.insert(1,pSW)
    #add midpoint to make sure some points are inflow labelled
    #inflow is on bottom
    #vertices.insert(2,[pSW[0]+0.5*inflowLength,pSW[1]])
    #vertices.insert(3,[pSW[0]+inflowLength,pSW[1]])
    #
    #vertices.append([bathymetryPoints[-1][0]+outflowPad-0.5*outflowLength,bathymetryPoints[-1][1]])
    #right
    #vertices.append(pSE)
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
    #segments.append([1,2])
    #segmentFlags.append(segmentLabels['inflow'])
    #segments.append([2,3])
    #segmentFlags.append(segmentLabels['inflow'])
    for i in range(1,nvertices-1):
        segments.append([i,i+1])
        segmentFlags.append(segmentLabels['bottom'])
    #segments.append([nvertices-4,nvertices-3])
    #segmentFlags.append(segmentLabels['outflow'])
    #segments.append([nvertices-3,nvertices-2])
    #segmentFlags.append(segmentLabels['outflow'])
    segments.append([nvertices-2,nvertices-1])
    segmentFlags.append(segmentLabels['right'])
    segments.append([nvertices-1,0])
    segmentFlags.append(segmentLabels['top'])
    print vertices,"s",segments,"sF",segmentFlags
    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags)
    domain.backHeight = backHeight
    domain.x_0 = x_0;  domain.x_be= x_be
    domain.x_bs= x_bs; 
    domain.x_ce= x_ce; domain.x_bse=x_bse
    domain.inflowLength= inflowLength; domain.beachLength = beachLength
    domain.beachSlope = beachSlope; domain.B = B; domain.backSlopeLength = backSlopeLength
    domain.backSlope  = backSlope
    domain.bathymetry = bathymetry ; domain.bathymetryGrad = bathymetryGrad
    domain.x_L=x_L ; domain.z_L = z_L; domain.L = L
    #go ahead and add a boundary tags member 
    domain.boundaryTags = segmentLabels
    return domain

if __name__=='__main__':
    import os
    domain =  lwi_dike_2d()
    domain.writeAsymptote("lwi_dike_2d")
    domain.writePoly("lwi_dike_2d")
    os.system("asy -V lwi_dike_2d")
