#! /usr/bin/env python
import proteus
from proteus.MeshTools import *
## \ingroup scripts
#
# \file genMesh.py
#
# \brief A script for generating (multilevel) meshes on cuboidal domains.

def genMesh():
    from optparse import OptionParser
    usage = "usage: %genMesh [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--boundaryMesh",
                      help="print boundary mesh to 3dm file",
                      action="store_true",
                      dest="boundaryMesh",
                      default=False)
    parser.add_option("-f", "--file",
                      help="write adh mesh to MESHFILE",
                      action="store",
                      type="string",
                      dest="meshFile",
                      default="adh")
    parser.add_option("-F", "--format",
                      help="write mesh in FORMAT",
                      action="store",
                      type="string",
                      dest="format",
                      default="adh")
    parser.add_option("-x", "--nnodes-x",
                      help="use NX nodes in the x direction",
                      action="store",
                      type="int",
                      dest="nx",
                      default=2)
    parser.add_option("-y", "--nnodes-y",
                      help="use NY nodes in the y direction",
                      action="store",
                      type="int",
                      dest="ny",
                      default=2)
    parser.add_option("-z", "--nnodes-z",
                      help="use NZ nodes in the z direction",
                      action="store",
                      type="int",
                      dest="nz",
                      default=2)
    parser.add_option("-X", "--length-x",
                      help="use LX units as the x dimension of the domain",
                      action="store",
                      type="float",
                      dest="lx",
                      default=1.0)
    parser.add_option("-Y", "--length-y",
                      help="use LY units as the y dimension of the domain",
                      action="store",
                      type="float",
                      dest="ly",
                      default=1.0)
    parser.add_option("-Z", "--length-z",
                      help="use LZ units as the z dimension of the domain",
                      action="store",
                      type="float",
                      dest="lz",
                      default=1.0)
    parser.add_option("-m", "--matlab",
                      help="print edges to files for plotting in matlab",
                      action="store_true",
                      dest="matlab",
                      default=False)
    parser.add_option("-M", "--view-matlab",
                      help="plot edges in matlab",
                      action="store_true",
                      dest="viewMatlab",
                      default=False)
    parser.add_option("-g", "--gnuplot",
                      help="print edges to files for plotting in gnuplot",
                      action="store_true",
                      dest="gnuplot",
                      default=False)
    parser.add_option("-G", "--view-gnuplot",
                      help="plot edges in gnuplot",
                      action="store_true",
                      dest="viewGnuplot",
                      default=False)
    parser.add_option("-r", "--refine",
                      help="refine mesh NR times",
                      action="store",
                      dest="nr",
                      type="int",
                      default=0)
    parser.add_option("-R", "--refine-method",
                      help="refine mesh using RTYPE refinement, RTYPE = 4T or FB",
                      action="store",
                      dest="rType",
                      type="string",
                      default='FB')
    parser.add_option("-c", "--cube-splitting-method",
                      help="split cubes using CTYPE refinement, CTYPE = 5T or 6T",
                      action="store",
                      dest="cType",
                      type="string",
                      default='6T')
    (opts, args) = parser.parse_args()
    print("The domain is Lx=%5.5e by Ly=%5.5e by Lz=%5.5e" % \
          (opts.lx,opts.ly,opts.lz))
    print("The mesh is nx=%i by ny=%i by nz=%i" % \
          (opts.nx,opts.ny,opts.nz))
    grid=RectangularGrid(opts.nx,opts.ny,opts.nz,opts.lx,opts.ly,opts.lz)
    cmesh=[]
    mesh=[]
    if opts.nz > 1:
        cmesh=TetrahedralMesh()
        if opts.cType == '5T':
            cmesh.rectangularToTetrahedral5T(grid)
        elif opts.cType == '6T':
            cmesh.rectangularToTetrahedral6T(grid)
        else:
            cmesh.rectangularToTetrahedral6T(grid)
        for i in range(opts.nr):
            mesh=TetrahedralMesh()
            if opts.rType == 'FB':
                mesh.refineFreudenthalBey(cmesh)
            elif opts.rType == '4T':
                mesh.refine4T(cmesh)
            else:
                print("refinement type not recognized, using 4T refinment")
                mesh.refine4T(cmesh)
            cmesh = mesh
            print(mesh.meshInfo())
        else:
            mesh = cmesh
    elif opts.ny > 1:
        print("Building triangular mesh")
        cmesh=TriangularMesh()
        cmesh.rectangularToTriangular(grid)
        for i in range(opts.nr):
            mesh=TriangularMesh()
            if opts.rType == 'FB':
                mesh.refineFreudenthalBey(cmesh)
            elif opts.rType == '4T':
                mesh.refine3t(cmesh)
            else:
                print("refinement type not recognized, using 4T refinment")
                mesh.refine3t(cmesh)
            cmesh = mesh
        else:
            mesh = cmesh
    else:
        cmesh=EdgeMesh()
        cmesh.rectangularToEdge(grid)
        for i in range(opts.nr):
            mesh=EdgeMesh()
            mesh.refine2e(cmesh)
            cmesh = mesh
        else:
            mesh = cmesh
    print("Writing ADH mesh to "+opts.meshFile)
    if opts.format == 'adh':
        mesh.writeMeshADH(opts.meshFile)
        if opts.boundaryMesh == True:
            mesh.boundaryMesh.writeMeshADH(opts.meshFile+'_boundary.3dm')
            mesh.writeBoundaryFacesADH('boundaryFaces_'+opts.meshFile+'.bc')
            mesh.writeBoundaryNodesADH('boundaryNodes_'+opts.meshFile+'.bc')
    elif opts.format == 'ensight':
        mesh.writeMeshEnsight(opts.meshFile)
    print(mesh.meshInfo())
    if opts.matlab or opts.viewMatlab:
        mesh.writeEdgesMatlab('meshEdges')
    if opts.viewMatlab:
        mesh.viewMeshMatlab('meshEdges')
    if opts.gnuplot or opts.viewGnuplot:
        mesh.writeEdgesGnuplot('meshEdges')
        mesh.boundaryMesh.writeEdgesGnuplot('boundaryMeshEdges')
    if opts.viewGnuplot:
        mesh.viewMeshGnuplot('meshEdges')

if __name__ == '__main__':
    #import profile
    #profile.run('genMesh()','genMeshProf')
    #import pstats
    #p = pstats.Stats('genMeshProf')
    #p.sort_stats('cumulative').print_stats(20)
    #p.sort_stats('time').print_stats(20)
    genMesh()
