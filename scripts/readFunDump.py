#! /usr/bin/env python
from numpy import *
import proteus
from proteus.MeshTools import *
from proteus.FemTools import *
## \ingroup scripts
#
# \file readFunDump.py
#
# \brief A script for reading an ADH .3dm mesh and .data file to build a finite element function.

def readFunLite():
    from optparse import OptionParser
    usage = "usage: %readFunLite [options] meshFile funFile funOutFile"
    parser = OptionParser(usage=usage)
 #     parser.add_option("-m", "--matlab",
#                        help="print edges to files for plotting in matlab",
#                        action="store_true",
#                        dest="matlab",
#                        default=False)
#      parser.add_option("-M", "--view-matlab",
#                        help="plot edges in matlab",
#                        action="store_true",
#                        dest="viewMatlab",
#                        default=False)
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
    (opts, args) = parser.parse_args()
    if len(args) == 4:
        meshFilename = args[1]
        funFilename = args[2]
        funOutFilename = args[3]
    else:
        print(usage)
        exit(1)
    #just get the nodes from the mesh file
    meshIn= open(meshFilename,'r')
    print("Reading nodes of the mesh from file=%s" % (meshFilename))
    line = meshIn.readline()
    columns = line.split()
    nodes=[]
    adhBase = 1
    while (columns[0] != 'ND'):
        line = meshIn.readline()
        columns = line.split()
    else:
        while (line and columns[0] == 'ND'):
            nodes.append(Node(int(columns[1]) - adhBase,
                              float(columns[2]),
                              float(columns[3]),
                              float(columns[4])))
            line=meshIn.readline()
            columns=line.split()
    funIn = open(funFilename,'r')
    nNodes=0
    fun = []
    for i in range(6):
        line = funIn.readline()
        print(line.strip())
        words = line.split()
        if words[0] == 'ND':
            nNodes = int(words[1])
            if nNodes != len(nodes):
                print("the number of nodes in mesh and function files don't match")
    line = funIn.readline()
    while line.strip() != 'ENDDS':
        print("Reading "+line.strip())
        words = line.split()
        u = zeros(nNodes,Float)
        for i in range(nNodes):
            u[i] = float(funIn.readline())
        fun.append(u)
        line = funIn.readline()
    print("Read %i timesteps" % len(fun))
    print("Writing coordinate files x.grf, y.grf, and z.grf")
    xiList={}
    xOut = open('x.grf','w')
    yOut = open('y.grf','w')
    zOut = open('z.grf','w')
    for n in nodes:
        xOut.write('%15.8e \n'% n.x)
        yOut.write('%15.8e \n'% n.y)
        zOut.write('%15.8e \n'% n.z)
    xOut.close()
    yOut.close()
    zOut.close()
    for ts,f in enumerate(fun):
        funOutFileTS = funOutFilename+repr(ts)+'.grf'
        print("Writing time step=%i to file=%s" % (ts,funOutFileTS))
        funOut = open(funOutFileTS,'w')
        for v in f:
            funOut.write('%15.8e \n'% v)
        funOut.close()
if __name__ == '__main__':
    import sys
    #import profile
    #profile.run('readFunLite(sys.argv[1],sys.argv[2],sys.argv[3])','readFunLiteProf')
    #import pstats
    #p = pstats.Stats('readFunLiteProf')
    #p.sort_stats('cumulative').print_stats(20)
    #p.sort_stats('time').print_stats(20)
    readFunLite()
