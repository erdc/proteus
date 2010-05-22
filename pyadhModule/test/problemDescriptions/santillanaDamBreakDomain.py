#! /usr/bin/env python
import math
from pyadh import Domain

#simple domain (rotated clockwise 90 degrees)
#   upstream
#    width
#   -------  
#   |     |  up_len
#   |     |  
#   --   --  
#     | |    mid_len
#     | |    
#   --   --  
#   |     |  down_len
#   |     |
#   -------
#   downstream

def twoReservoirDomain(up_width=3.0,
                       mid_width=1.0,
                       up_len=2.0,
                       mid_len=2.0,
                       down_len=2.0,
                       name="santillanaDamBreakDomain"):
    #label walls separately for bcs
    boundaries=['upstream','downstream','wall_horiz','wall_vert']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    vertexBoundaries = ['upstream','downstream','wall']
    vertexBoundaryTags=dict([(key,i+1) for (i,key) in enumerate(vertexBoundaries)])


    flume_len = up_len+mid_len+down_len
    L = (flume_len,up_width)
    mid_y_start=(up_width-mid_width)*0.5
    mid_y_end  = (up_width+mid_width)*0.5
    segments = []; segmentFlags = []; vertices = []; vertexFlags = []
    #work aroung going counter clockwise
    vertices.append([0.0,0.0]); vertexFlags.append(vertexBoundaryTags['upstream'])
    vertices += [[up_len,0.0],
                [up_len,mid_y_start],
                [up_len+mid_len,mid_y_start],
                [up_len+mid_len,0.0]]
    for i in range(1,len(vertices)):
        vertexFlags.append(vertexBoundaryTags['wall'])
    vertices.append([flume_len,0.0]); vertexFlags.append(vertexBoundaryTags['downstream'])

    nv = len(vertices)
    for i in range(nv-1):
        segments.append([i,i+1]); 

    segmentFlags.extend([boundaryTags['wall_horiz'],
                         boundaryTags['wall_vert'],
                         boundaryTags['wall_horiz'],
                         boundaryTags['wall_vert'],
                         boundaryTags['wall_horiz']])
                        
    #end bottom boundary
    #now downstream
    vertices.append([flume_len,up_width]); vertexFlags.append(vertexBoundaryTags['downstream'])
    segments.append([nv-1,nv]); segmentFlags.append(boundaryTags['downstream'])
    #now top
    vertices += [[flume_len-down_len,up_width],
                [flume_len-down_len,mid_y_end],
                [up_len,mid_y_end],
                [up_len,up_width]]
    vertexFlags.extend([vertexBoundaryTags['wall'],
                        vertexBoundaryTags['wall'],
                        vertexBoundaryTags['wall'],
                        vertexBoundaryTags['wall']])
        
    for i in range(nv,len(vertices)):
        segments.append([i,i+1]);
    segmentFlags.extend([boundaryTags['wall_horiz'],
                         boundaryTags['wall_vert'],
                         boundaryTags['wall_horiz'],
                         boundaryTags['wall_vert'],
                         boundaryTags['wall_horiz']])
    #finish with upstream
    vertices.append([0.0,up_width]); vertexFlags.append(vertexBoundaryTags['upstream'])
    segments.append([len(vertices)-1,0]); segmentFlags.append(boundaryTags['upstream'])

    domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                  vertexFlags=vertexFlags,
                                                  segments=segments,
                                                  segmentFlags=segmentFlags,
                                                  name=name)

    domain.boundaryTags = boundaryTags
    return domain
#     name='dummy'
#     poly = open(name+'.poly','w')
#     poly.write('%d %d %d %d \n' % (len(vertices),2,0,1))
#     poly.write("#vertices \n")
#     for vN,v in enumerate(vertices):
#         poly.write('%d %12.5e %12.5e %d #%s \n' % (vN+1,v[1][0],v[1][1],v[2],v[0]))
#     poly.write('%d %d \n' % (len(vertices),1))
#     poly.write("#segments \n")
#     for sN,s in enumerate(segments):
#         poly.write('%d %d %d %d \n' % (sN+1,s[0]+1,s[1]+1,s[2]))
#     poly.write("#segments \n")
#     poly.write('%d\n' % (0,))
#     poly.write("#regions \n")
#     poly.write('%d\n' % (1,))
#     poly.write('%d %12.5e %12.5e %d' % (1,
#                                         vertices[0][1][0]+up_len*1.0e-8,
#                                         vertices[0][1][1]+up_len*1.0e-8,
#                                         0+1))
#     poly.close()
#     return boundaryTags,L
if __name__=='__main__':
    import os
    domain = twoReservoirDomain(name="santillanaDamBreakDomain")
    domain.writeAsymptote("santillanaDamBreakDomain")
    domain.writePoly("santillanaDamBreakDomain")
    os.system("asy -V santillanaDamBreakDomain")
    os.system('triangle -Aq30den santillanaDamBreakDomain.poly')
    os.system('showme santillanaDamBreakDomain.1.ele')
