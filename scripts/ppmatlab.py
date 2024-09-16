## Automatically adapted for numpy.oldnumeric Apr 14, 2008 by -c

def writeMeshMatlabFormat(mesh,meshFileBase):
    """
    build array data structures for matlab finite element mesh representation
    and write to a file to view and play with in matlatb

    in matlab can then print mesh with

    pdemesh(p,e,t)

    where

    p is the vertex or point matrix
    e is the edge matrix, and
    t is the element matrix

    points matrix is [2 x num vertices]
    format :
    row 1 = x coord,
    row 2 = y coord for nodes in mesh

    edge matrix is [7 x num edges]
    format:
    row 1 = start vertex number
    row 2 = end vertex number
    row 3 = start value in edge parameterization, should be 0
    row 4 = end   value in edge parameterization, should be 1
    row 5 = global edge id, base 1
    row 6 = subdomain on left? always 1 for now
    row 7 = subdomain on right? always 0 for now

    element matrix is [4 x num elements]
    row 1 = vertex 1 global number
    row 2 = vertex 2 global number
    row 3 = vertex 3 global number
    row 4 = triangle subdomain number
    where 1,2,3 is a local counter clockwise numbering of vertices in
    triangle

    """
    import numpy as numpy
    matlabBase = 1
    p = numpy.zeros((2,mesh['nNodes_global']),numpy.float_)
    e = numpy.zeros((7,mesh['nElementBoundaries_global']),numpy.float_)
    t = numpy.zeros((4,mesh['nElements_global']),numpy.float_)

    #load p,e,t and write file
    mfile = open(meshFileBase+'.m','w')

    mfile.write('p = [ ... \n')
    for nN in range(mesh['nNodes_global']):
        p[0,nN]=mesh['nodeArray'][nN,0]
        p[1,nN]=mesh['nodeArray'][nN,1]
        mfile.write('%g %g \n' % tuple(p[:,nN]))
    mfile.write(']; \n')
    mfile.write("p = p\';\n")  #need transpose for matlab

    mfile.write('e = [ ... \n')
    for ebN in range(mesh['nElementBoundaries_global']):
        e[0,ebN]=mesh['elementBoundaryNodesArray'][ebN,0] + matlabBase #global node number of start node base 1
        e[1,ebN]=mesh['elementBoundaryNodesArray'][ebN,1] + matlabBase #global node number of end node base 1
        e[2,ebN]=0.0 #edge param. is 0 to 1
        e[3,ebN]=1.0
        e[4,ebN]=ebN + matlabBase  #global edge number base 1
        e[5,ebN]=0 #subdomain to left
        e[6,ebN]=1 #subdomain to right
        mfile.write('%g %g %g %g %g %g %g \n' % tuple(e[:,ebN]))
    mfile.write(']; \n')
    mfile.write("e = e\';\n")  #need transpose for matlab

    #write triangles last
    mfile.write('t = [ ... \n')
    for eN in range(mesh['nElements_global']):
        t[0,eN]=mesh['elementNodesArray'][eN,0]+matlabBase    #global node number for vertex 0
        t[1,eN]=mesh['elementNodesArray'][eN,1]+matlabBase    #global node number for vertex 0
        t[2,eN]=mesh['elementNodesArray'][eN,2]+matlabBase    #global node number for vertex 0
        t[3,eN]=1                     #subdomain id
        mfile.write('%g %g %g %g \n' % tuple(t[:,eN]))
    mfile.write(']; \n');
    mfile.write("t = t\';\n") #need transpose for matlab


    mfile.close()
    return p,e,t


########################################################################
if __name__ == '__main__':
    import os,shelve
    import ppmatlab,numpy.oldnumeric as numpy

    os.listdir('./results')

    filename = './results/re_forsyth2_ss_2d_pre_forsyth2_ss_2d_c0p1_n_mesh_results.dat'

    res = shelve.open(filename)

    mesh = res['mesh']

    mmfile = 'forsyth2MeshMatlab'
    p,e,t = ppmatlab.writeMeshMatlabFormat(mesh,mmfile)
