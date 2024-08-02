import math
n_domain_vertices = 100
domain_vertices =[(0.75*math.sin(2.0*math.pi*float(n)/float(n_domain_vertices))+0.5,0.75*math.cos(2.0*math.pi*float(n)/float(n_domain_vertices))+0.5) for  n in range(n_domain_vertices)]
grain_centers  = [(0.25,0.25),(0.25,0.75),(0.75,0.75),(0.75,0.25),(0.5,0.5),(0.5+0.5,0.5),(0.5-0.5,0.5),(0.5,0.5-0.5),(0.5,0.5+0.5)]
radii = [0.15,0.15,0.15,0.15,0.125,0.1,0.1,0.1,0.1]
points_on_grain = 50
nvertices = len(domain_vertices) + len(grain_centers)*points_on_grain
poly = open('pome.poly','w')
poly.write('%d %d %d %d \n' % (nvertices,2,0,1))
#write vertices
poly.write("#vertices \n")
for v,p in enumerate(domain_vertices):
    poly.write('%d %18.12e %18.12e %d \n' % (v+1,p[0],p[1],1))
#write segments
nSegmentGroups = 0
segments=[]
segmentGroups=[]
for sN in range(len(domain_vertices)-1):
    segments.append([sN,sN+1])
    segmentGroups.append(nSegmentGroups)
segments.append([len(domain_vertices)-1,0])
segmentGroups.append(nSegmentGroups)
vStart = len(domain_vertices)
sStart = len(segments)
nSegmentGroups = nSegmentGroups+1
for g,c in enumerate(grain_centers):
    for gb in range(points_on_grain):
        pb  = (radii[g]*math.sin(float(gb)/float(points_on_grain)*2.0*math.pi),radii[g]*math.cos(float(gb)/float(points_on_grain)*2.0*math.pi))
        poly.write('%d %18.12e %18.12e %d \n' % (vStart + gb+1,c[0]+pb[0],c[1]+pb[1],1+g))
    for gb in range(points_on_grain-1):
        segments.append([sStart+gb,sStart+gb+1])
        segmentGroups.append(nSegmentGroups)
    segments.append([sStart+points_on_grain-1,sStart])
    segmentGroups.append(nSegmentGroups)
    vStart = vStart + points_on_grain
    sStart = sStart + points_on_grain
    nSegmentGroups = nSegmentGroups+1
poly.write('%d %d \n' % (len(segments),1))
print(segments)
poly.write("#segments \n")
for sN,s,sG in zip(list(range(len(segments))),segments,segmentGroups):
    poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,sG))
poly.write('%d \n' % len(grain_centers))
for gN,g in enumerate(grain_centers):
    poly.write('%d %18.12e %18.12e \n' % (gN,g[0],g[1]))
poly.close()
