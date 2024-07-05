#!/usr/bin/env python
nBlocks = 22
blockLeft  = [0.0, 0.7/3.0, 0.0, 0.0, 0.05, 0.1, 0.2, 0.5, 0.6, 0.2, 0.45, 0.65, 0.25, 0.1, 0.05, 0.1, 0.5, 0.1, 0.2, 0.35, 0.1, 0.35]
blockRight = [0.7, 2.0*0.7/3.0,0.7, 0.05, 0.1, 0.2, 0.45, 0.65, 0.65, 0.5, 0.5, 0.7, 0.35, 0.6, 0.65, 0.2, 0.6, 0.25, 0.5, 0.6, 0.6, 0.6]
blockFront = [0.0, 0.45, 0.0, 0.05, 0.05, 0.05, 0.1, 0.05, 0.15, 0.05, 0.1, 0.05, 0.2, 0.3, 0.4, 0.35, 0.35, 0.2, 0.35, 0.2, 0.15, 0.25]
blockBack  = [0.5, 0.5, 0.05, 0.5 , 0.4, 0.15, 0.15, 0.15, 0.4, 0.1, 0.15, 0.5, 0.3, 0.35, 0.5, 0.4, 0.4, 0.3, 0.4, 0.25, 0.2, 0.3]
vertices = set()
for bN in range(nBlocks):
    vertices |= set([(blockLeft[bN],blockFront[bN]),
                     (blockLeft[bN],blockBack[bN]),
                     (blockRight[bN],blockFront[bN]),
                     (blockRight[bN],blockBack[bN])])
nVertices = len(vertices)
c2v = dict([(xy,vN) for vN,xy in enumerate(vertices)])
poly = open('KueperFrind.poly','w')
poly.write('#PSLG representation of Kueper and Frind problem \n')
poly.write('%d %d %d %d \n' % (nVertices,2,0,1))
poly.write('#vertices \n')
for vN,xy in enumerate(vertices):
    poly.write('%d %12.5e %12.5e %d\n' % (vN+1,xy[0],xy[1],1))
poly.write('%d %d\n' % (nBlocks*4,1))
poly.write('#segments\n')
for bN in range(nBlocks):
    poly.write('%d %d %d %d\n' % (1+bN*4+ 0,
                                  1+c2v[(blockLeft[bN],blockFront[bN])],
                                  1+c2v[(blockRight[bN],blockFront[bN])],
                                  1+bN))
    poly.write('%d %d %d %d\n' % (1+bN*4+ 1,
                                  1+c2v[(blockRight[bN],blockFront[bN])],
                                  1+c2v[(blockRight[bN],blockBack[bN])],
                                  1+bN))
    poly.write('%d %d %d %d\n' % (1+bN*4+ 2,
                                  1+c2v[(blockRight[bN],blockBack[bN])],
                                  1+c2v[(blockLeft[bN],blockBack[bN])],
                                  1+bN))
    poly.write('%d %d %d %d\n' % (1+bN*4+ 3,
                                  1+c2v[(blockLeft[bN],blockBack[bN])],
                                  1+c2v[(blockLeft[bN],blockFront[bN])],
                                  1+bN))
poly.write('#holes\n0\n#regions\n0\n')
poly.close()
