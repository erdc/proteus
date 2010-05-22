#!/usr/bin/env python

def genPoly(fileprefix,
            depth=0.7/10.0,
            ny=2,
            source_x = (0.7/3.0,0.7*2.0/3.0),
            source_y = (0.0,0.7/10.0)):
    nBlocks = 22
    blockLeft  = [0.00, 0.00, 0.05, 0.10, 0.20, 0.50, 0.60, 0.20, 0.45, 0.65, 0.25, 0.10, 0.05, 0.10, 0.50, 0.10, 0.20, 0.35, 0.10, 0.35]
    blockRight = [0.70, 0.05, 0.10, 0.20, 0.45, 0.65, 0.65, 0.50, 0.50, 0.70, 0.35, 0.60, 0.65, 0.20, 0.60, 0.25, 0.50, 0.60, 0.60, 0.60]
    blockBottom = [0.00, 0.05, 0.05, 0.05, 0.10, 0.05, 0.15, 0.05, 0.10, 0.05, 0.20, 0.30, 0.40, 0.35, 0.35, 0.20, 0.35, 0.20, 0.15, 0.25]
    blockTop  = [0.05, 0.50, 0.40, 0.15, 0.15, 0.15, 0.40, 0.10, 0.15, 0.50, 0.30, 0.35, 0.50, 0.40, 0.40, 0.30, 0.40, 0.25, 0.20, 0.30]
#     blockLeft  = [0.0, 0.7/3.0, 0.0, 0.0, 0.05, 0.1, 0.2, 0.5, 0.6, 0.2, 
#                   0.45, 0.65, 0.25, 0.1, 0.05, 0.1, 0.5, 0.1, 0.2, 0.35, 0.1, 0.35]
#     blockRight = [0.7, 2.0*0.7/3.0,0.7, 0.05, 0.1, 0.2, 0.45, 0.65, 0.65, 
#                   0.5, 0.5, 0.7, 0.35, 0.6, 0.65, 0.2, 0.6, 0.25, 0.5, 0.6, 0.6, 0.6]
#     blockBottom = [0.0, 0.45, 0.0, 0.05, 0.05, 0.05, 0.1, 0.05, 0.15, 0.05, 
#                   0.1, 0.05, 0.2, 0.3, 0.4, 0.35, 0.35, 0.2, 0.35, 0.2, 0.15, 0.25]
#     blockTop  = [0.5, 0.5, 0.05, 0.5 , 0.4, 0.15, 0.15, 0.15, 0.4, 0.1, 
#                   0.15, 0.5, 0.3, 0.35, 0.5, 0.4, 0.4, 0.3, 0.4, 0.25, 0.2, 0.3]
    xList = set(blockLeft)
    xList |= set(blockRight)
    xList = list(xList)
    xList.sort()
    #
    zList = set(blockBottom)
    zList |= set(blockTop)
    zList = list(zList)
    zList.sort()
    dy = depth/float(ny-1)
    yList = [i*dy for i in range(ny)]
    boundaries=['left',
                'right',
                'front',
                'back',
                'bottom',
                'top',
                'source']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])
    nx = len(xList)
    ny = len(yList)
    nz = len(zList)
    nxy = nx*ny
    nxyz = nxy*nz
    poly = open(fileprefix+'.poly','w')
    poly.write('#PSLG representation of Kueper and Frind problem, extruded to 3d \n')
    poly.write('%d %d %d %d \n' % (nxyz,3,0,1))
    poly.write('#vertices \n')
    def nN(i,j,k):
        return i*nxy+j*nx+k+1
    def bt(i,j,k):
        if i==0:
            return boundaryTags['bottom']
        if i==nz-1:
            (x,y) = (xList[k],yList[j])
            if source_x[0] <= x  and x <= source_x[1] and source_y[0] <= y and y <= source_y[1]:
                return boundaryTags['source']
            else:
                return boundaryTags['top']
        if j==0:
            return boundaryTags['front']
        if j==ny-1:
            return boundaryTags['back']
        if k==0:
            return boundaryTags['left']
        if k==nx-1:
            return boundaryTags['right']
        return 0
    def bt_topFacets(i,j,k):
        if i==0:
            return boundaryTags['bottom']
        if i==nz-1:
            isSource = True
            for jj in range(2):
                for kk in range(2):
                    (x,y) = (xList[k+kk],yList[j+jj])
                    if not (source_x[0] <= x  and x <= source_x[1] and source_y[0] <= y and y <= source_y[1]):
                        isSource=False
            if isSource:
                return boundaryTags['source']
            else:
                return boundaryTags['top']
        return 0
    for i,z in enumerate(zList):
        for j,y in enumerate(yList):
            for k,x in enumerate(xList):
                poly.write('%d %12.5e %12.5e %12.5e %d\n' % (nN(i,j,k),x,y,z,bt(i,j,k)))
    poly.write('#facets\n')
    ncells = nx*(ny-1)*(nz-1)+ny*(nx-1)*(nz-1)+nz*(nx-1)*(ny-1)
    poly.write('%d %d\n' % (ncells,1))
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                if j < ny-1 and k < nx-1:
                    poly.write('1 0 %d\n' % (bt_topFacets(i,j,k),))
                    poly.write('4 %d %d %d %d\n' % (nN(i,j,k),
                                                    nN(i,j+1,k),
                                                    nN(i,j+1,k+1),
                                                    nN(i,j,k+1)))
                if i < nz-1 and k < nx-1:
                    poly.write('1 0 %d\n' % (bt(-1,j,-1),))
                    poly.write('4 %d %d %d %d\n' % (nN(i,j,k),
                                                    nN(i+1,j,k),
                                                    nN(i+1,j,k+1),
                                                    nN(i,j,k+1)))
                if j < ny-1 and i < nz-1:
                    poly.write('1 0 %d\n' % (bt(-1,-1,k),))
                    poly.write('4 %d %d %d %d\n' % (nN(i,j,k),
                                                    nN(i+1,j,k),
                                                    nN(i+1,j+1,k),
                                                    nN(i,j+1,k)))
    poly.write('#holes\n0\n')
    poly.write('#regions\n%d\n' % ((nx-1)*(nz-1),))
    yc = 0.5*depth
    rN=1
    for i in range(nz-1):
        for k in range(nx-1):
            (xc,zc) = (0.5*(xList[k]+xList[k+1]),0.5*(zList[i]+zList[i+1]))
            for nB in range(len(blockLeft)):
                if (xc > blockLeft[nB] and
                    xc < blockRight[nB] and
                    zc > blockBottom[nB] and
                    zc < blockTop[nB]):
                    break
            poly.write('%d %12.5e %12.5e %12.5e %d\n' % (rN,xc,yc,zc,nB))
            rN+=1
    poly.close()
    return boundaryTags

if __name__=='__main__':
    import os
    top=0.5
    right=0.7
    genPoly('kueper3d_test',depth=0.7,ny=4,source_y = (0.7/3.0,0.7*2.0/3.0))
    #triangleOptions = "pqAfena%e" % ((top*right*top*0.1/6.0)*0.001,)
    #print triangleOptions
    #os.system('tetgen -'+triangleOptions+' kueper3d_test.poly')
    #os.system('tetview kueper3d_test.1.ele')
    from pyadh import Domain
    d = Domain.PiecewiseLinearComplexDomain("kueper3d_test")
    d.writeAsymptote("kueper3d_test.asy")
    os.system("asy -V kueper3d_test.asy")
