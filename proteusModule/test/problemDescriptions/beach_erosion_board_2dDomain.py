#! /usr/bin/env python
import math
from pyadh import Domain

"""
 generate waves and runup on  WES beach erosion board embankment

 units should be meters by default


 z=0         beach slope
             (1/10)         1/s slope           1/s_2 slope
 |--------|----------------|-----------|-------|-----------|-----|
 0        x_bs            x_be         x_se   x_ce        x_bse  L

x_bs = inflow pad (2 [m] default), 
z(x) = 0, [0,x_bs]

x_be = x_bs + beachLength (4.48 [m])
z(x) = 0.0 + (x-x_bs)beachSlope, [x_bs,x_be]

x_se = x_be + (h_c + h_s)*s, h_c and h_s from Sitanggang Lynnett report
z(x) = (x-x_be)1/s + z(x_be), [x_be,x_se] 

x_ce = x_se + B [m]
z(x) = z(x_se), [x_se,x_ce]

x_bse= x_ce + (h_c + h_s)*0.5*s_2
z(x) = z(x_ce) + (x-x_ce)*1/s_2 [x_ce,x_bse]

x_L  = x_bse + outflowLength  [m]
z(x) = z(x_bse), [x_bse,L]

total domain height = z(x_se) + domainHeightPad
"""

def beach_erosion_board_2d(domainHeightPad=2.0,
                           inflowLength=2.0,
                           beachLength=4.48,
                           beachSlope =0.1,
                           h_c=0.054,
                           h_s=0.081,
                           s  =3.0,
                           B  = 0.0896,
                           s_2= -0.2,
                           backStepFraction=0.25,
                           outflowLength=0.5):

    """
     generate waves and runup on  WES beach erosion board embankment

     units should be meters by default


     z=0         beach slope
                 (1/10)         1/s slope           1/s_2 slope
     |--------|----------------|-----------|-------|-----------|-----|
     0        x_bs            x_be         x_se   x_ce        x_bse  L

    x_bs = inflow pad (2 [m] default), 
    z(x) = 0, [0,x_bs]

    x_be = x_bs + beachLength (4.48 [m])
    z(x) = 0.0 + (x-x_bs)beachSlope, [x_bs,x_be]

    x_se = x_be + (h_c + h_s)*s, h_c and h_s from Sitanggang Lynnett report
    z(x) = (x-x_be)1/s + z(x_be), [x_be,x_se] 

    x_ce = x_se + B [m]
    z(x) = z(x_se), [x_se,x_ce]

    x_bse= x_ce + (h_c + h_s)*backStepFraction*s_2
    z(x) = z(x_ce) + (x-x_ce)*1/s_2 [x_ce,x_bse]

    x_L  = x_bse + outflowLength  [m]
    z(x) = z(x_bse), [x_bse,L]

    total domain height = z(x_se) + domainHeightPad
    """
    #describe bathymetry as above
    x_0  = 0.0
    x_bs = x_0 + inflowLength; x_be = x_bs + beachLength
    x_se = x_be + (h_c+h_s)*s
    x_ce = x_se + B
    x_bse= x_ce + (h_c+h_s)*backStepFraction*abs(s_2)

    def bathymetry(x):
        if x[0] <= x_bs:  return 0.0
        if x[0] <= x_be:  return 0.0 + (x[0]-x_bs)*beachSlope
        if x[0] <= x_se:  return bathymetry([x_be]) + (x[0]-x_be)/s
        if x[0] <= x_ce:  return bathymetry([x_se])
        if x[0] <= x_bse: return bathymetry([x_ce]) + (x[0]-x_ce)/s_2
        return bathymetry([x_bse])
    def bathymetryGrad(x):
        if x[0] <= x_bs:  return (0.0,)
        if x[0] <= x_be:  return (beachSlope,) #beach slope
        if x[0] <= x_se:  return (1.0/s,)
        if x[0] <= x_ce:  return (0.0,)
        if x[0] <= x_bse: return (1./s_2,)
        return (0.0,)
    
    x_L  = x_bse + outflowLength
    z_L  = bathymetry([x_se]) + domainHeightPad
    L = [x_L,z_L]

    xpoints = [x_0,x_bs,x_be,x_se,x_ce,x_bse,x_L]
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
    domain.x_bs= x_bs; domain.x_se= x_se
    domain.x_ce= x_ce; domain.x_bse=x_bse
    domain.inflowLength= inflowLength; domain.beachLength = beachLength
    domain.s = s; domain.B = B; domain.s_2 = s_2
    domain.h_c= h_c; domain.h_s = h_s
    domain.bathymetry = bathymetry ; domain.bathymetryGrad = bathymetryGrad
    #go ahead and add a boundary tags member 
    domain.boundaryTags = segmentLabels
    return domain

if __name__=='__main__':
    import os
    domain =  beach_erosion_board_2d()
    domain.writeAsymptote("beach_erosion_board2d")
    domain.writePoly("beach_erosion_board2d")
    os.system("asy -V beach_erosion_board2d")
