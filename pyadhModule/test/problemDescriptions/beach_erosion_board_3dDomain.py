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

def beach_erosion_board_3d(domainHeightPad=2.0,
                           inflowLength=2.0,
                           beachLength=4.48,
                           beachSlope =0.1,
                           h_c=0.054,
                           h_s=0.081,
                           s  =3.0,
                           B  = 0.0896,
                           s_2= -0.2,
                           backStepFraction=0.25,
                           outflowLength=0.5,
                           width = 0.5):

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
        if x[0] <= x_bs:  return (0.0,0.0)
        if x[0] <= x_be:  return (beachSlope,0.0) #beach slope
        if x[0] <= x_se:  return (1.0/s,0.0)
        if x[0] <= x_ce:  return (0.0,0.0)
        if x[0] <= x_bse: return (1./s_2,0.0)
        return (0.0,0.0)
    
    x_L  = x_bse + outflowLength
    z_L  = bathymetry([x_se]) + domainHeightPad
    L = [x_L,width,z_L]

    xpoints = [x_0,x_bs,x_be,x_se,x_ce,x_bse,x_L]
    zpoints = [bathymetry([x]) for x in xpoints]
    bathymetryPoints = [[x,0.0,z] for x,z in zip(xpoints,zpoints)]

    backHeight = bathymetry([x_bse])

    #have to assume points in correct order

    #pad for inflow outflow
    pSW = bathymetryPoints[0]
    pSE= bathymetryPoints[-1]

    #now get top corners of domain
    minY = 0.0
    pNW = [pSW[0],0.0,minY+z_L]
    pNE = [pSE[0],0.0,minY+z_L]

    
    #check vertical coordinates
    tmp = sorted(bathymetryPoints,cmp=lambda x,y: int(x[2]-y[2]))    
    assert minY <= tmp[0][1], "found point below proposed block floor minY=%s tmp[0]= " % (minY,tmp[0])
    assert minY+ z_L > tmp[-1][1], "found point above proposed block ceiling maxnY=%s tmp[-1]= " % (minY+z_L,
                                                                                                    tmp[-1])
    frontVertices = [p for p in bathymetryPoints]

    #start with NW corner and work way around
    #left 
    frontVertices.insert(0,pNW)

    frontVertices.append(pNE)
    nfrontVertices = len(frontVertices)

    backVertices = [[v[0],width,v[2]] for v in frontVertices]
    vertices = frontVertices+backVertices


    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

    facets = []
    facetFlags=[]
    front_facets=[[vN for vN in range(nfrontVertices)]]
    facets.append(front_facets)
 
    facetFlags.append(boundaryTags['front'])
    back_facets=[[vN for vN in range(nfrontVertices,2*nfrontVertices)]]
    facets.append(back_facets)
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
    for vN in range(1,nfrontVertices-2):
        bottom_facet=[[vN,vN+1,vN+1+nfrontVertices,vN+nfrontVertices]]
        facets.append(bottom_facet)
        facetFlags.append(boundaryTags['bottom'])
 
    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 facets=facets,
                                                 facetFlags=facetFlags)
    domain.backHeight = backHeight
 
    domain.backHeight = backHeight
    domain.x_0 = x_0;  domain.x_be= x_be
    domain.x_bs= x_bs; domain.x_se= x_se
    domain.x_ce= x_ce; domain.x_bse=x_bse
    domain.inflowLength= inflowLength; domain.beachLength = beachLength
    domain.s = s; domain.B = B; domain.s_2 = s_2
    domain.h_c= h_c; domain.h_s = h_s
    domain.bathymetry = bathymetry ; domain.bathymetryGrad = bathymetryGrad
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain


def beach_erosion_board_porous_3d(domainHeightPad=2.0,
                                  inflowLength=2.0,
                                  beachLength=4.48,
                                  beachSlope =0.1,
                                  h_c=0.054,
                                  h_s=0.081,
                                  s  =3.0,
                                  B  = 0.0896,
                                  s_2= -0.2,
                                  backStepFraction=0.25,
                                  outflowLength=0.5,
                                  width = 0.5):

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

    model the region between x_be and x_L as a homogeneous porous region
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
        if x[0] <= x_bs:  return (0.0,0.0)
        if x[0] <= x_be:  return (beachSlope,0.0) #beach slope
        if x[0] <= x_se:  return (1.0/s,0.0)
        if x[0] <= x_ce:  return (0.0,0.0)
        if x[0] <= x_bse: return (1./s_2,0.0)
        return (0.0,0.0)
    
    x_L  = x_bse + outflowLength
    z_L  = bathymetry([x_se]) + domainHeightPad
    L = [x_L,width,z_L]

    xpoints = [x_0,x_bs,x_be,x_se,x_ce,x_bse,x_L]
    zpoints = [bathymetry([x]) for x in xpoints]
    bathymetryPoints = [[x,0.0,z] for x,z in zip(xpoints,zpoints)]

    backHeight = bathymetry([x_bse])

    #have to assume points in correct order

    #pad for inflow outflow
    pSW = bathymetryPoints[0]
    pSE= bathymetryPoints[-1]

    #now get top corners of domain
    minY = 0.0
    pNW = [pSW[0],0.0,minY+z_L]
    pNE = [pSE[0],0.0,minY+z_L]


    #check vertical coordinates
    tmp = sorted(bathymetryPoints,cmp=lambda x,y: int(x[2]-y[2]))    
    assert minY <= tmp[0][1], "found point below proposed block floor minY=%s tmp[0]= " % (minY,tmp[0])
    assert minY+ z_L > tmp[-1][1], "found point above proposed block ceiling maxnY=%s tmp[-1]= " % (minY+z_L,
                                                                                                    tmp[-1])
    #add porous region
    xporousPoints = [x_be,x_se,x_ce,x_bse,x_L]
    zporousPoints = [bathymetry([x_be]) for x in xporousPoints]
    porousPoints=[[x,0.0,z] for x,z in zip(xporousPoints,zporousPoints)]

    frontVertices = [p for p in bathymetryPoints]

    #start with NW corner and work way around
    #left 
    frontVertices.insert(0,pNW)

    frontVertices.append(pNE)
    nfrontFluidVertices = len(frontVertices)

    frontVertices += [p for p in porousPoints]
    nfrontVertices = len(frontVertices)

    backVertices = [[v[0],width,v[2]] for v in frontVertices]
    vertices = frontVertices+backVertices


    boundaries=['left','right','bottom','top','front','back']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

    facets = []
    facetFlags=[]
    front_facets=[[vN for vN in range(nfrontFluidVertices)]]
    facets.append(front_facets)
    facetFlags.append(boundaryTags['front'])
    front_porous_facets = [[vN for vN in range(nfrontFluidVertices,nfrontVertices)]]
    for v in [nfrontFluidVertices-1,nfrontFluidVertices-2,nfrontFluidVertices-3,nfrontFluidVertices-4,nfrontFluidVertices-5]:
        front_porous_facets[-1].append(v)
    facets.append(front_porous_facets)
 
    facetFlags.append(boundaryTags['front'])
    back_facets=[[vN for vN in range(nfrontVertices,nfrontVertices+nfrontFluidVertices)]]
    facets.append(back_facets)
    facetFlags.append(boundaryTags['back'])
    back_porous_facets=[[vN for vN in range(nfrontVertices+nfrontFluidVertices,2*nfrontVertices)]]
    top_facets=[[vN for vN in [0,nfrontFluidVertices-1,nfrontVertices+nfrontFluidVertices-1,nfrontVertices]]]#[[vN for vN in range(0,nfrontVertices+1)]]
    facets.append(top_facets)
    facetFlags.append(boundaryTags['top'])
    right_facets=[[vN for vN in [nfrontFluidVertices-2,nfrontFluidVertices-1,nfrontVertices+nfrontFluidVertices-1,nfrontVertices+nfrontFluidVertices-2]]]#[[vN for vN in range(1,1+(nfrontVertices+1))]]
    right_facets+=[[vN for vN in [nfrontVertices-1,nfrontFluidVertices-2,nfrontVertices+nfrontFluidVertices-2,2*nfrontVertices-1]]]#[[vN for vN in range(1,1+(nfrontVertices+1))]]
    facets.append(right_facets)
    facetFlags.append(boundaryTags['right'])
    left_facets=[[vN for vN in [0,nfrontVertices,nfrontVertices+1,1]]]
    facets.append(left_facets)
    facetFlags.append(boundaryTags['left'])
    for vN in range(1,nfrontFluidVertices-2):
        bottom_facet=[[vN,vN+1,vN+1+nfrontVertices,vN+nfrontVertices]]
        facets.append(bottom_facet)
        facetFlags.append(boundaryTags['bottom'])
    #
    #add a region for porous domain
    xporousCenter = sum(xporousPoints)/len(xporousPoints)
    zporousCenter = 0.5*(bathymetry([x_be])+bathymetry([x_L]))
    regions = [[xporousCenter,0.5*width,zporousCenter]]
    regionFlags = [1]

    domain = Domain.PiecewiseLinearComplexDomain(vertices=vertices,
                                                 facets=facets,
                                                 facetFlags=facetFlags,
                                                 regions=regions,
                                                 regionFlags=regionFlags)
    domain.backHeight = backHeight
 
    domain.backHeight = backHeight
    domain.x_0 = x_0;  domain.x_be= x_be
    domain.x_bs= x_bs; domain.x_se= x_se
    domain.x_ce= x_ce; domain.x_bse=x_bse
    domain.inflowLength= inflowLength; domain.beachLength = beachLength
    domain.s = s; domain.B = B; domain.s_2 = s_2
    domain.h_c= h_c; domain.h_s = h_s
    domain.bathymetry = bathymetry ; domain.bathymetryGrad = bathymetryGrad
    #go ahead and add a boundary tags member 
    domain.boundaryTags = boundaryTags
    return domain

if __name__=='__main__':
    import os
    porousDomain = True
    if porousDomain == False:
        domain =  beach_erosion_board_3d()
        domain.writeAsymptote("beach_erosion_board_3d")
        domain.writePoly("beach_erosion_board_3d")
        os.system("asy -V -render=5 beach_erosion_board_3d")
        os.system("asy -f png -render 3 beach_erosion_board_3d")
        os.system("asy -prc -f pdf -render=5 beach_erosion_board_3d")
    else:
        domain =  beach_erosion_board_porous_3d()
        domain.writeAsymptote("beach_erosion_board_porous_3d")
        domain.writePoly("beach_erosion_board_porous_3d")
        os.system("asy -V -render=5 beach_erosion_board_porous_3d")
        os.system("asy -f png -render 3 beach_erosion_board_porous_3d")
        os.system("asy -prc -f pdf -render=5 beach_erosion_board_porous_3d")
