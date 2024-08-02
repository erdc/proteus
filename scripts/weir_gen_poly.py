import math
vertices = [ ('bottom_left' , (0.0,0.0)),
             ('weir_bottom_left' , (0.8 + 1.5 + 6.66,0.0)),
             ('weir_top_left' , (0.8 + 1.5 + 6.66,0.235)),
             ('weir_top_right' , (0.8 + 1.5 + 6.66+0.02,0.235)),
             ('weir_bottom_right' , (0.8 + 1.5 + 6.66+0.02,0.0)),
             ('bottom_right' , (0.8 + 1.5 + 6.66 + 5.7 + 6.74,0.0)),
             ('top_right', (0.8 + 1.5 + 6.66 + 5.7 + 6.74,0.61)),
             ('top_left' , (0.0,0.61))]
# vertices = [ ('bottom_left' , (0.0,0.0)),
#              ('weir_bottom_left' , (0.8 + 1.5 + 6.66,0.0)),
#              ('weir_top_left' , (0.8 + 1.5 + 6.66,0.235)),
#              ('weir_top_right' , (0.8 + 1.5 + 6.66+1.02,0.235)),
#              ('weir_bottom_right' , (0.8 + 1.5 + 6.66+1.02,0.0)),
#              ('bottom_right' , (0.8 + 1.5 + 6.66 + 5.7 + 6.74,0.0)),
#              ('top_right', (0.8 + 1.5 + 6.66 + 5.7 + 6.74,0.61)),
#              ('top_left' , (0.0,0.61))]
nvertices = len(vertices)
poly = open('weir.poly','w')
poly.write('%d %d %d %d \n' % (nvertices,2,0,1))
#write vertices
poly.write("#vertices \n")
for v,p in enumerate(vertices):
    poly.write('%d %12.5e %12.5e %d #%s \n' % (v+1,p[1][0],p[1][1],1,p[0]))
#write segments
nSegmentGroups = 0
segments=[]
segmentGroups=[]
for sN in range(nvertices-1):
    segments.append([sN,sN+1])
    segmentGroups.append(nSegmentGroups)
segments.append([nvertices-1,0])
segmentGroups.append(nSegmentGroups)
poly.write('%d %d \n' % (nvertices,1))
poly.write("#segments \n")
for sN,s,sG in zip(list(range(len(segments))),segments,segmentGroups):
    poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,sG))
poly.write('%d' % (0,))
poly.close()
