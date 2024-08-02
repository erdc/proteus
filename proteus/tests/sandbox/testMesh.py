#!/usr/bin/env python
from MeshTools import *
import numpy

def constructTriangularMeshOnRectangle(Lx,Ly,nx,ny,writeMesh=0,
                                       meshFileBase='mesh2d'):
    """
    wrapper function for making a triangular mesh on the rectangle
    [0,Lx] x [0,Ly].

    viewMesh is a flag to allow printing mesh when constructed
    viewMesh -- 0 no visualization
                1 gnuplot
                2 matlab
    """
    nz = 1
    Lz = 1.0
    grid2d = RectangularGrid(nx,ny,nz,Lx,Ly,Lz)
    #grid2d.writeEdgesGnuplot('grid2d')
    #grid2d.viewMeshGnuplotPipe('grid2d')

    mesh2d = TriangularMesh()

    mesh2d.rectangularToTriangular(grid2d)

    if writeMesh == 1:
        #print mesh in gnuplot format
        mesh2d.writeEdgesGnuplot(meshFileBase)
        #can view with
        #mesh2d.viewMeshGnuplotPipe(meshFileBase)
    elif writeMesh == 2:
        mesh2d.writeEdgesMatlab(meshFileBase)
        #view in matlab with meshFileBase.m
    #end else

    return mesh2d
#end construct triangular mesh

def buildMatlabMeshDataStructures(mesh,meshFileBase='meshMatlab'):
    """
    build array data structures for matlab finite element mesh representation
    and write to a file to view and play with in matlatb

    in matlab can then print mesh with

      pdemesh(p,e,t)

    where p is the vertex or point matrix
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
    matlabBase = 1
    #total number of nodes
    nnodes = mesh.nNodes_global
    mnodes = Numeric.zeros((2,nnodes),Numeric.Float)
    for key in list(mesh.nodeDict.keys()):
        nN = mesh.nodeDict[key].N      #global node number?
        xy = mesh.nodeDict[key].p      #physical location
        mnodes[0,nN] = xy[X]
        mnodes[1,nN] = xy[Y]
    #end mnodes build loop

    #now get edge array
    nedges = mesh.nEdges_global
    #mwf debug
    #print 'in build edge number of edges= ',nedges
    #print 'in build edge number of nodes= ',nnodes
    #print 'in build edge number of elems= ',mesh.nElements_global
    medges = Numeric.zeros((7,nedges),Numeric.Int)

    for key in list(mesh.edgeDict.keys()):
        nE = mesh.edgeDict[key].N     #global edge number
        pE = mesh.edgeDict[key].nodes #nodes defining edge
        #mwf debug
        #print 'in build edge ',nE,' nodes ',pE[0],' ',pE[0].N,' ',pE[1],' ',pE[1].N
        medges[0,nE] = pE[0].N+matlabBase        #global node number
                                                 #of start node base 1
        medges[1,nE] = pE[1].N+matlabBase        #global node number
                                                 #of end node base 1
        medges[2,nE] = 0.0            #edge param. is 0 to 1
        medges[3,nE] = 1.0
        medges[4,nE] = nE+matlabBase  #global edge number base 1
        medges[5,nE] = 0              #subdomain to left?
        medges[6,nE] = 1              #subdomain to right?
    #end

    #now get element array
    nelems = mesh.nElements_global
    melems = Numeric.zeros((4,nelems),Numeric.Int)

    for key in list(mesh.triangleDict.keys()):
        nE = mesh.triangleDict[key].N        #global element number
        pE = mesh.triangleDict[key].nodes    #nodes defining triangle

        #global numbers here need to be base 1
        melems[0,nE] = pE[0].N+matlabBase    #global node number for vertex 0
        melems[1,nE] = pE[1].N+matlabBase    #global node number for vertex 1
        melems[2,nE] = pE[2].N+matlabBase    #global node number for vertex 2
        melems[3,nE] = 1                     #subdomain id
    #end elements loop

    #now write to a file
    mfile = open(meshFileBase+'.m','w')

    mfile.write('p = [ ... \n')
    for j in range(mnodes.shape[1]):
        sn = '%g %g \n' % (mnodes[0,j],mnodes[1,j])
        mfile.write(sn)
    #end nodes loop
    mfile.write(']; \n')
    mfile.write("p = p\';\n")  #need transpose for matlab

    mfile.write('e = [ ... \n')
    for j in range(medges.shape[1]):
        sn = '%g %g %g %g %g %g %g \n' % tuple(medges[:,j])
        mfile.write(sn)
    #end edge loop
    mfile.write(']; \n')
    mfile.write("e = e\';\n")  #need transpose for matlab

    #write triangles last
    mfile.write('t = [ ... \n')
    for j in range(melems.shape[1]):
        sn = '%g %g %g %g \n' % tuple(melems[:,j])
        mfile.write(sn)
    #end element loop
    mfile.write(']; \n');
    mfile.write("t = t\';\n") #need transpose for matlab


    mfile.close()
    return (mnodes,medges,melems)

#end build matlab data structures



if __name__=='__main__':

    #how much junk to print out
    verboseLevel = 2
    #first just create a simple triangular mesh and look at it in a
    #couple of different ways
    Lx = 1.0   #domain length in x and y
    Ly = 1.0

    #number of nodes for rectangular grid upon which triangular mesh
    #will be built should get 2 triangles for each rectangle
    #(nx-1)(ny-1) in the original grid
    nx = 3
    ny = 3

    #flag for viewing mesh in construction
    #0 -- do nothing (default)
    #1 -- gnuplot
    #2 -- matlab
    viewMesh = 2
    mesh = constructTriangularMeshOnRectangle(Lx,Ly,nx,ny,viewMesh)

    print('mesh Info says \n',mesh.meshInfo())
    fileName2 = 'meshV2'
    mp,me,mt = buildMatlabMeshDataStructures(mesh,fileName2)

    if verboseLevel > 1:
        #do brute force loop through array to look at it
        print('matlab node array is ')
        for j in range(mp.shape[1]): #number of columns is number of nodes
            print('\t',mp[0,j],' ',mp[1,j])
        #end for

        #do brute force loop through edge array too
        print('matlab edge array holds (matlab edge id, node 0, node 1)')
        print('note base 0')
        for j in range(me.shape[1]): #number of columns is number of edges
            print('\t',me[4,j]-1,' ',me[0,j]-1,' ',me[1,j]-1)
        #end for

        #do brute force loop through element array too
        print('matlab elem array (matlab elem id, node 0, node 1, node 3)')
        print('note base 0')
        for j in range(mt.shape[1]): #number of columns is number of edges
            print('\t',j,' ',mt[0,j]-1,' ',mt[1,j]-1,' ',mt[2,j]-1)
        #end for
    #end verbose print out for mesh
