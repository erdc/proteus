"""
A class hierarchy and tools for building domains of PDE's.

.. inheritance-diagram:: proteus.Domain
   :parts: 1
"""

import numpy as np
from proteus.MeshTools import MeshParallelPartitioningTypes as mpt


class D_base:
    """
    The base class for domains
    """
    def __init__(self, nd, name="defaultDomain", units="m"):
        """
        Set dimensions (nd), name string, and units string
        """
        if nd not in [1,2,3]:
            raise RuntimeError("Domain object must have dimension 1,2, or 3")
        self.nd = nd
        self.name = name
        self.units = units
        # make default empty list for parameters
        self.vertices = []
        self.vertexFlags = []
        self.segments = []
        self.segmentFlags = []
        self.facets = []  # necessary in 2D for gmsh
        self.facetFlags = []
        self.holes = []
        self.holes_ind = []  # for gmsh: index of hole (2D: facet; 3D: volume)
        self.regions = []
        self.regionFlags = []
        self.volumes = []  # only for gmsh: loop of facets
        # for bounding box
        self.x = []
        self.L = []
        # list of shape instances (from proteus.SpatialTools.py)
        self.shape_list = []
        # list of boundary conditions
        self.bc = []
        # list of auxiliaryVariables
        self.auxiliaryVariables = []
        # needed for gmsh
        self.volumes = []
        self.boundaryTags = {}
        # attach a MeshOption class
        self.MeshOptions = MeshOptions(self.nd)

    def writeAsymptote(self, fileprefix):
        """
        Write a representation of the domain to file using the Asymptote vector graphics language.
        """
        raise UserWarning("Asymptote output not implemented")
    def writeXdmf(self,ar):
        """
        Store the domain in an Xdmf file.
        """
        raise UserWarning("Xdmf output not implemented")

    def getBoundingBox(self):
        """
        Get the bounding box of a set of vertices.
        """
        if self.x or self.L:
            self.x = []
            self.L = []
        for d in range(self.nd):
            xMax = max(row[d] for row in self.vertices)
            xMin = min(row[d] for row in self.vertices)
            self.x += [xMin]
            self.L += [xMax-xMin]


    def writeGeo(self, fileprefix, group_names=False):
        """
        Write the PSLG using gmsh geo format incomplete write now
        probably run into problems with orientation, no concept of a
        line loop dummyAxis (0,1,2) is the remaining axis for
        embedding the 2d geometry in 3d and xref is the constant value
        for that axis
        """
        self.geofile = fileprefix+'.geo'
        self.polyfile = fileprefix
        geo = open(self.geofile,'w')
        pp = {}  # physical points
        pl = {}  # physical lines
        ps = {}  # physical surfaces
        pv = {}  # physical volume
        sN = len(self.segments)

        # Vertices
        print self.vertices
        geo.write('\n// Points\n')
        z = 0
        for i, v in enumerate(self.vertices):
            if self.nd == 3:
                z = v[2]
            geo.write("Point(%d) = {%g,%g,%g};\n" % (i+1,v[0],v[1], z))
            if self.vertexFlags:
                flag = self.vertexFlags[i]
                if flag in pp:
                    pp[flag] += [i+1]
                else:
                    pp[flag] = [i+1]
        nb_points = i+1

        lines_dict = {}
        for i in range(nb_points):
            lines_dict[i] = {}

        # Lines
        geo.write('\n// Lines\n')
        # line_list = []
        for i, s in enumerate(self.segments):
            geo.write("Line(%d) = {%d,%d};\n" % (i+1,s[0]+1,s[1]+1))
            # add segments in dictionary
            lines_dict[s[0]][s[1]] = i
            if self.segmentFlags:
                flag = self.segmentFlags[i]
                if flag in pl:
                    pl[flag] += [i+1]
                else:
                    pl[flag] = [i+1]

        # Surfaces
        geo.write('\n// Surfaces\n')
        lines = 0
        lineloop_count = 0
        surface_line_list = []  # need to store for mesh constraints later
        # facet
        for i, f in enumerate(self.facets):
            seg_flag = sN+i+1
            lineloop = []
            lineloops = []
            lineloops_list = []
            line_list = []
            # subfacet
            if self.nd == 3 or (self.nd == 2 and i not in self.holes_ind):
                for j, subf in enumerate(f):
                    lineloop = []
                    # vertices in facet
                    print 'subf'+str(subf)
                    for k, ver in enumerate(subf):
                        if ver in lines_dict[subf[k-1]].keys():
                            lineloop += [lines_dict[subf[k-1]][ver]+1]
                        elif subf[k-1] in lines_dict[ver].keys():
                            # reversed
                            lineloop += [-(lines_dict[ver][subf[k-1]]+1)]
                        else:
                            ind = seg_flag+lines
                            lines += 1
                            geo.write('Line(%d) = {%d,%d};\n' % (ind, subf[k-1]+1, ver+1))
                            lineloop += [ind]
                    line_list += lineloop
                    geo.write('Line Loop(%d) = {%s};\n' % (lineloop_count+1, str(lineloop)[1:-1]))
                    lineloops += [lineloop_count+1]
                    lineloop_count += 1
                    # print i, j
                surface_line_list += [line_list]
                geo.write('Plane Surface(%d) = {%s};\n' % (i+1, str(lineloops)[1:-1]))
                if self.facetFlags:
                    flag = self.facetFlags[i]
                    if flag in ps:
                        ps[flag] += [i+1]
                    else:
                        ps[flag] = [i+1]

        # Volumes
        geo.write('\n// Volumes\n')
        for i, V in enumerate(self.volumes):
            surface_loops = []
            if i not in self.holes_ind:
                for j, sV in enumerate(V):
                    lineloop_count += 1
                    geo.write('Surface Loop(%d) = {%s};\n' % (lineloop_count, str((np.array(sV)+1).tolist())[1:-1]))
                    surface_loops += [lineloop_count]
                geo.write('Volume(%d) = {%s};\n' % (i+1, str(surface_loops)[1:-1]))
                flag = self.regionFlags[i]
                if self.regionFlags:
                    if flag in pv:
                        pv[flag] += [i+1]
                    else:
                        pv[flag] = [i+1]

        # Physical Groups
        geo.write('\n// Physical Groups\n')
        if self.boundaryTags:
            inv_bt = {v: k for k, v in self.boundaryTags.iteritems()}
        for flag in pp:
            ind = pp[flag]
            if self.boundaryTags and group_names is True:
                flag = '"'+inv_bt[flag]+'", '+str(flag)
            geo.write("Physical Point({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))
        for flag in pl:
            ind = pl[flag]
            if self.boundaryTags and group_names is True:
                flag = '"'+inv_bt[flag]+'", '+str(flag)
            geo.write("Physical Line({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))
        for flag in ps:
            ind = ps[flag]
            if self.boundaryTags and group_names is True:
                flag = '"'+inv_bt[flag]+'", '+str(flag)
            geo.write("Physical Surface({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))
        for flag in pv:
            ind = pv[flag]
            if self.boundaryTags and group_names is True:
                flag = '"'+inv_bt[flag]+'", '+str(flag)
            geo.write("Physical Volume({0}) = {{{1}}};\n".format(str(flag), str(ind)[1:-1]))



class MeshOptions:
    """
    Mesh options for the domain
    """
    def __init__(self, nd):
        self.nd = nd
        self.he = None
        self.genMesh = True
        self.use_gmsh = False
        self.outputFiles = {'name': 'mesh',
                            'poly': True,
                            'ply': False,
                            'asymptote': False}
        self.restrictFineSolutionToAllMeshes = False
        self.parallelPartitioningType = mpt.node
        self.nLayersOfOverlapForParallel = 0
        if self.nd == 2:
            self.triangle_string = 'VApq30Dena'
        if self.nd == 3:
            self.triangle_string = 'VApq1.35q12feena'
        self.triangleOptions = None  # defined when setTriangleOptions called

    def elementSize(self, he, refinement_lvl=0.):
        self.he = he*0.5**refinement_lvl

    def parallelPartitioningType(self, type='node', layers_overlap=0):
        if type == 'element' or 0:
            self.parallelPartitioningType = mpt.element
        if type == 'node' or 1:
            self.parallelPartitioningType = mpt.node
        self.nLayersOfOverlapForParallel = layers_overlap

    def setTriangleOptions(self):
        if self.he is None:
            self.he = 1.
        if self.nd == 2:
            self.triangleOptions = self.triangle_string + '%8.8f' \
                                   % (self.he**2/2.,)
        elif self.nd == 3:
            self.triangleOptions = self.triangle_string + '%21.16e' \
                                   % (self.he**3/6.,)
    def outputFiles(self, name='mesh', poly=True, ply=False, asymptote=False):
        self.outputFiles['name'] = name
        self.outputFiles['poly'] = poly
        self.outputFiles['ply'] = ply
        self.outputFiles['asymptote'] = asymptote


class RectangularDomain(D_base):
    """
    Rectangular domains in 1,2, or 3D specified by the lower left hand corner and dimensions.
    """
    def __init__(self,
                 L=[1.0,1.0,1.0],
                 x=[0.0,0.0,0.0],
                 name="DefaultRectangularDomain",
                 units="m"):
        D_base.__init__(self,len(L),name,units)
        self.polyfile = None
        self.boundaryTags = {'left':3,
                             'right':5,
                             'front':2,
                             'back':6,
                             'top':4,
                             'bottom':1}
        if self.nd==2:
            self.boundaryTags = {'left':4,
                                 'right':2,
                                 'top':3,
                                 'bottom':1}

        self.x=x
        self.L=L
    def writePoly(self,fileprefix):
        """
        Write the RectangularDomain using the poly format.
        """
        self.polyfile = fileprefix #[temp] see PSLG for better implementation that checks if it already exists
        self.boundaryLegend = self.boundaryTags
        unitesize=4.0/self.L[0]
        f = open(fileprefix+".poly",'w')
        if self.nd==2:
            self.boundaryLegend = self.boundaryTags
            fileString="""
# vertices
4 2 0 1
1 %(x0)f %(x1)f 1
2 %(x0pL0)f %(x1)f 1
3 %(x0pL0)f %(x1pL1)f 3
4 %(x0)f %(x1pL1)f 3
# segments
4 1
1 1 2 1
2 2 3 2
3 3 4 3
4 4 1 4
# holes
0
# regions
1
1 %(x0pL0)f %(x1pL1)f 1 """  % {'x0': self.x[0],'x1': self.x[1],'x0pL0': self.x[0]+self.L[0],'x1pL1': self.x[1]+self.L[1]}
# tjp altered the regions line and the boundary tags from 0 to 1 to get to work
        elif self.nd==3:
            fileString="""
# vertices
8 3 0 0
1 %(x0)f %(x1)f %(x2)f 1
2 %(x0)f %(x1pL1)f %(x2)f 2
3 %(x0pL0)f %(x1pL1)f %(x2)f 3
4 %(x0pL0)f %(x1)f %(x2)f 4
5 %(x0)f %(x1)f %(x2pL2)f 5
6 %(x0)f %(x1pL1)f %(x2pL2)f 6
7 %(x0pL0)f %(x1pL1)f %(x2pL2)f 7
8 %(x0pL0)f %(x1)f %(x2pL2)f 8
# facets
6 1
1 0 1
4 1 2 3 4                       # bottom
1 0 2
4 1 5 8 4                       # front
1 0 3
4 1 5 6 2                       # left
1 0 4
4 5 6 7 8                       # top
1 0 5
4 4 8 7 3                       # right
1 0 6
4 3 7 6 2                       # back
# holes
0
# regions
1
1 %(m0)f %(m1)f %(m2)f 1 """ % {'x0': self.x[0],
                                'x1': self.x[1],
                                'x2': self.x[2],
                                'x0pL0': self.x[0]+self.L[0],
                                'x1pL1': self.x[1]+self.L[1],
                                'x2pL2': self.x[2]+self.L[2],
                                'm0':(self.x[0]+self.L[0])/2.0,
                                'm1':(self.x[1]+self.L[1])/2.0,
                                'm2':(self.x[2]+self.L[2])/2.0}
        f.write(fileString)
        f.close()
    def writeAsymptote(self, fileprefix):
        """
        Write the rectangular domain as a box displaying the dimesions and coloring the boundaries according to the boundary flags.
        """
        unitsize=4.0/self.L[0]
        f = open(fileprefix+".asy",'w')
        if self.nd==1:
            fileString="""
unitsize(4.0 inches/%(L)f);
size(5 inches);
real L=%(L)f;
real offset=.0125L;
real x=%(x)f;
string str="$%(L)2.2f\mbox{%(units)s}$";
import graph;
import palette;
pen[] allPens=Wheel();
pen[] myPens = new pen[3];
for(int i=0;i< 3;++i)
  {
   int iPen = round(i*allPens.length/3);
   myPens[i] = allPens[iPen];
  }
draw((x,0)--(x+L,0)^^(x,offset)--(x+L,offset),myPens[0]);
draw((x,0)--(x,offset),myPens[1]);
draw((x+L,0)--(x+L,offset),myPens[2]);
draw(str,(x,-offset)--(x+L,-offset),black,Bars,Arrows,PenMargins);
""" % {'L':self.L[0],'x':self.x[0],'units':self.units}
        elif self.nd==2:
            fileString="""unitsize(4.0 inches / %(Lx)f);
size(5 inches);
real Lx=%(Lx)f;
real Ly=%(Ly)f;
real offset=0.0125Lx;
real x=%(x)f;
real y=%(y)f;
string strx="$%(Lx)2.2f\mbox{%(units)s}$";
string stry="$%(Ly)2.2f\mbox{%(units)s}$";
import graph;
import palette;
pen[] allPens=Wheel();
pen[] myPens = new pen[4];
for(int i=0;i< 4;++i)
  {
   int iPen = round(i*allPens.length/4);
   myPens[i] = allPens[iPen];
  }
draw((x,y)--(x+Lx,y),myPens[0]);
draw((x,y+Ly)--(x+Lx,y+Ly),myPens[1]);
draw((x,y)--(x,y+Ly),myPens[2]);
draw((x+Lx,y)--(x+Lx,y+Ly),myPens[3]);
draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);
draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);
""" % {'Lx':self.L[0],'Ly':self.L[1],'x':self.x[0],'y':self.x[1],'units':self.units}
        elif self.nd==3:
            fileString="""import three;
currentprojection=perspective(-2,-2,1,up=Z,target=O,showtarget=true,autoadjust=true,center=false);
unitsize(4.0 inches/%(Lx)f);
size(5 inches);
real Lx=%(Lx)f;
real Ly=%(Ly)f;
real Lz=%(Lz)f;
real offset=.0125Lx;
real x=%(x)f;
real y=%(y)f;
real z=%(z)f;
triple p0=(x,y,z);
triple p1=(x,y+Ly,z);
triple p2=(x+Lx,y+Ly,z);
triple p3=(x+Lx,y,z);
triple p4=(x,y,z+Lz);
triple p5=(x,y+Ly,z+Lz);
triple p6=(x+Lx,y+Ly,z+Lz);
triple p7=(x+Lx,y,z+Lz);
string strx="$%(Lx)2.2f\mbox{%(units)s}$";
string stry="$%(Ly)2.2f\mbox{%(units)s}$";
string strz="$%(Lz)2.2f\mbox{%(units)s}$";
import graph;
import palette;
pen[] allPens=Wheel();
pen[] myPens = new pen[6];
for(int i=0;i< 6;++i)
  {
   int iPen = round(i*allPens.length/6);
   myPens[i] = allPens[iPen];
  }
draw(surface(p0--p1--p2--p3--cycle),myPens[0]);
draw(surface(p4--p5--p6--p7--cycle),myPens[1]);
draw(surface(p0--p1--p5--p4--cycle),myPens[2]);
draw(surface(p1--p2--p6--p5--cycle),myPens[3]);
draw(surface(p2--p3--p7--p6--cycle),myPens[4]);
draw(surface(p3--p0--p4--p7--cycle),myPens[5]);
draw(strx,(x,y-offset,z-offset)--(x+Lx,y-offset,z-offset),-(Y+Z),black,Bars3(Z),Arrows3);
draw(stry,(x-offset,y,z-offset)--(x-offset,y+Ly,z-offset),-(X+Z),black,Bars3(Z),Arrows3);
draw(strz,(x-offset,y-offset,z)--(x-offset,y-offset,z+Lz),-(X+Y),black,Bars3(X),Arrows3);
shipout();
""" % {'Lx':self.L[0],'Ly':self.L[1],'Lz':self.L[2],'x':self.x[0],'y':self.x[1],'z':self.x[2],'units':self.units}
        f.write(fileString)
        f.close()
    def writeXdmf(self,ar):
        raise UserWarning("Xdmf output not implemented")

class GMSH_3D_Domain(D_base):
    """
    A domain described by a gmsh .geo file

    Parameters
    ----------
    geofile : str
              The base filename of the main .geo file
    he : float
         The maximum element diameter (this should go way)
    name : str
           A name for the domain
    units : str
            The units of length
    """
    def __init__(self, geofile, he,
                 name="GMSH_Domain",
                 units="m",
                 length_scale=1.0,
                 permute_dims=[0,1,2]):
        D_base.__init__(self, 2, name, units)
        self.geofile=geofile+".geo"
        self.polyfile=geofile
        self.he = he
        self.length_scale=length_scale
        self.permute_dims=permute_dims

class PlanarStraightLineGraphDomain(D_base):
    """
    2D domains described by planar straight line graphs.
    """

    def __init__(self, fileprefix=None, vertices=None, segments=None, holes=None, regions=None, vertexFlags=None,
                 segmentFlags=None, regionFlags=None, regionConstraints=None, bc=None, name="DefaultPSLGDomain", units="m"):
        """
        Construct the PSLG from lists of vertices, segments, etc. If no vertex or segment flags are given, then they are assigned as zero.
        """
        D_base.__init__(self, 2, name, units)
        if fileprefix:
            self.readPoly(fileprefix)
        else:
            self.polyfile = None
            self.vertices = vertices or []
            self.segments = segments or []
            self.holes = holes or []
            self.regions = regions or []
            self.vertexFlags = vertexFlags or []
            self.segmentFlags = segmentFlags or []
            self.regionFlags = regionFlags or []
            self.regionConstraints = regionConstraints or []
            self.update()

    def update(self):
            if self.regionConstraints:
                assert (self.regionFlags)
            if self.vertices:
                self.getBoundingBox()
            if self.segmentFlags:
                self.getSegmentPartition()
            else:
                if self.segments:
                    self.sMin = 0
                    self.sMax = len(self.segments) - 1
                    self.segmentFlags = [i for i in range(self.sMax + 1)]

    def getSegmentPartition(self):
        self.sMin=self.sMax=self.segmentFlags[0]
        for s in self.segmentFlags:
            if s > self.sMax:
                self.sMax=s
            if s < self.sMin:
                self.sMin=s

    def readPoly(self,fileprefix):
        """
        Read in the PSLG using the poly format.
        """
        self.polyfile = fileprefix
        self.name=fileprefix
        f = open(fileprefix+".poly",'r')
        firstLine=f.readline().split()
        while len(firstLine) == 0 or firstLine[0][0] == '#':
            firstLine = f.readline().split()
        nVertices = int(firstLine[0])
        self.dim = int(firstLine[1])
        assert(self.dim == 2)
        nVertexAttributes = int(firstLine[2])
        self.vertexAttributes = [[] for j in range(nVertexAttributes)]
        hasVertexFlag = int(firstLine[3])
        self.vertices=[]
        self.base=None
        self.vertexFlags = None
        if hasVertexFlag:
            self.vertexFlags=[]
        for i in range(nVertices):
            line = f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            if not self.base:
                self.base = int(line[0])
            self.vertices.append([float(line[1]),float(line[2])])
            for j in range(nVertexAttributes):
                self.vertexAttributes[j].append(float(line[3+j]))
            if hasVertexFlag:
                self.vertexFlags.append(line[3+nVertexAttributes])
        segmentLine = f.readline().split()
        while len(segmentLine) == 0 or segmentLine[0][0] == '#':
            segmentLine = f.readline().split()
        nSegments = int(segmentLine[0])
        hasSegmentFlag = bool(int(segmentLine[1]))
        self.segments=[]
        self.segmentFlags = None
        if hasSegmentFlag:
            self.segmentFlags=[]
        for i in range(nSegments):
            line = f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            self.segments.append([int(line[1])-self.base,int(line[2])-self.base])
            if hasSegmentFlag:
                self.segmentFlags.append(int(line[3]))
        holeLine = f.readline().split()
        while len(holeLine) == 0 or holeLine[0][0] == '#':
            holeLine = f.readline().split()
        nHoles = int(holeLine[0])
        self.holes=[]
        for i in range(nHoles):
            line =  f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            self.holes.append([float(line[0]),float(line[1])])
        regionLine = f.readline().split()
        while len(regionLine) == 0 or regionLine[0][0] == '#':
            regionLine = f.readline().split()
        nRegions = int(regionLine[0])
        self.regions=[]
        self.regionFlags=[]
        self.areaConstraints=[]
        for i in range(nRegions):
            line =  f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            self.regions.append([float(line[1]),float(line[2])])
            if len(line) > 3:
                self.regionFlags.append(int(line[3]))
                if len(line) > 4:
                    if line[4][0] != '#':
                        self.areaConstraints.append(float(line[4]))
        self.getBoundingBox()
        if self.segmentFlags :
            self.getSegmentPartition()
        f.close()
    def writePoly(self,fileprefix):
        """
        Write the PSLG using the poly format.
        """
        if not self.polyfile:
            self.polyfile = fileprefix
            pf = open(fileprefix+'.poly','w')
            #write first line of poly file
            if self.vertexFlags :
                hasVertexFlags=1
            else:
                hasVertexFlags=0
            if self.segmentFlags :
                hasSegmentFlags=1
            else:
                hasSegmentFlags=0
            pf.write("%d %d %d %d \n" % (len(self.vertices),2,0,hasVertexFlags))
            #write the vertices
            for vN,v in enumerate(self.vertices):
                pf.write('%d %21.16e %21.16e ' % (vN+1,v[0],v[1]))
                #import pdb; pdb.set_trace()
                if self.vertexFlags :#write vertex flag if we have vertexFlags
                    pf.write('%d\n' % (self.vertexFlags[vN],))
                else:
                    pf.write('\n')
            #write the segments
            pf.write('%d %d \n' % (len(self.segments),hasSegmentFlags))
            for sN,s in enumerate(self.segments):
                pf.write('%d %d %d ' % (sN+1,s[0]+1,s[1]+1))
                if self.segmentFlags :
                    pf.write('%d \n' % (self.segmentFlags[sN],))
                else:
                    pf.write('\n')
            if self.holes :
                pf.write('%d\n' % (len(self.holes),))
                for hN,h in enumerate(self.holes):
                    pf.write('%d %21.16e %21.16e \n' % (hN+1,
                                                      h[0],
                                                      h[1]))
            else:
                pf.write('%d\n' % (0,))
            if self.regions :
                pf.write('%d\n' % (len(self.regions),))
                for rN,r in enumerate(self.regions):
                    pf.write('%d %21.16e %21.16e ' % (rN+1,
                                                   r[0],
                                                   r[1]))
                    if self.regionConstraints:
                        pf.write('%d %21.16e\n' % (self.regionFlags[rN],self.regionConstraints[rN]))
                    elif self.regionFlags:
                        pf.write('%d\n' % (self.regionFlags[rN]))
                    else:
                        pf.write('\n')
            else:
                pf.write('%d\n' % (0,))
            pf.close()
        else:
            print "File already exists, not writing polyfile: " +`self.polyfile`
    def writeAsymptote(self,fileprefix):
        """
        Write a representation of the PSLG in the Asymptote vector graphics language
        """
        unitsize=4.0/self.L[0]
        f = open(fileprefix+".asy",'w')
        fileString="""
unitsize(4.0 inches / %(Lx)f);
size(11 inches);
real Lx=%(Lx)f;
real Ly=%(Ly)f;
real offset=0.0125Lx;
real x=%(x)f;
real y=%(y)f;
string strx="$%(Lx)2.2f\mbox{%(units)s}$";
string stry="$%(Ly)2.2f\mbox{%(units)s}$";
draw(strx,(x,y-offset)--(x+Lx,y-offset),S,black,Bars,Arrows,PenMargins);
draw(stry,(x-offset,y)--(x-offset,y+Ly),W,black,Bars,Arrows,PenMargins);
import graph;
import palette;
pen[] allPens=Wheel();
pen[] myPens = new pen[%(nSegmentFlags)d+1];
for(int i=0;i<%(nSegmentFlags)d+1;++i)
  {
   int iPen = round(i*allPens.length/(%(nSegmentFlags)d+1));
   myPens[i] = allPens[iPen];
  }
""" % {'Lx':self.L[0],'Ly':self.L[1],'x':self.x[0],'y':self.x[1],'units':self.units,
       'nSegmentFlags':(self.sMax-self.sMin)}
        #now loop over segments
        for s,sFlag in zip(self.segments,self.segmentFlags):
            fileString+="draw((%f,%f)--(%f,%f),myPens[%d]+linewidth(0.01));\n" % tuple(self.vertices[s[0]]+self.vertices[s[1]]+[sFlag-self.sMin])
        #for vN,v in enumerate(self.vertices):
        #    fileString+="label(\"(%1.1f,%1.1f)\",(%f,%f),N);\n" % (v[0],v[1],v[0],v[1])
        f.write(fileString)
        f.close()
    def writePLY(self,fileprefix):
        """
        Write the PLC domain in the Stanford PLY format.
        """
        if True:#overwrite
            pf = open(fileprefix+'.ply','w')
            if self.vertexFlags :
                hasVertexFlags=1
            else:
                hasVertexFlags=0
            if self.segmentFlags :
                hasSegmentFlags=1
            else:
                hasSegmentFlags=0
            hasSegmentFlags=False
            hasVertexFlags=False
            pf.write("""ply
format ascii 1.0
comment author: Proteus
comment object: %s
""" % (self.name,))
            pf.write("""element vertex %d
property float x
property float y
property float z
""" % (len(self.vertices),))
            if hasVertexFlags:
                pf.write("property int flag\n")
            pf.write("element face %d \n" % (len(self.segments),))
            pf.write("property list uchar int vertex_index\n")
            if hasSegmentFlags:
                pf.write("property int flag\n")
            pf.write("end_header\n")
            #write the vertices
            for vN,v in enumerate(self.vertices):
                pf.write('%21.16e %21.16e %21.16e' % (v[0],v[1],0.0))
                if hasVertexFlags:
                    pf.write(" %d\n" % (self.vertexFlags[vN],))
                else:
                    pf.write("\n")
            #write the segments
            for sN,s in enumerate(self.segments):
                pf.write(`len(s)`)
                for vN in s:
                    pf.write(" "+`vN`)
                if hasSegmentFlags:
                    pf.write(" %d\n" % (self.segmentFlags[fN],))
                else:
                    pf.write("\n")
            pf.close()
        else:
            print "File already exists, not writing polyfile: " +`self.polyfile`
    def writeXdmf(self,ar):
        """
        Store the PSLG domain in an XDMF file. For now we store the information on holes in and Information element.
        """
        raise UserWarning("Xdmf output not implemented")

class InterpolatedBathymetryDomain(PlanarStraightLineGraphDomain):
    """
    2D domains with an attached bathymetry data set.
    Meshes are generated by locally refining the PSLG until the interpolation error for a bathymetry point cloud meets a prescribed tolerance.

    Intial support measures interpolation error in the :math:`L_{\infty}` norm.
    """
    def __init__(self,fileprefix=None,vertices=None,segments=None,
                 holes=None,regions=None,vertexFlags=None,segmentFlags=None,regionFlags=None,
                 name="DefaultInterpolatedBathymetry",units="m",
                 bathy=None,
                 bathyGridDim=None):
        PlanarStraightLineGraphDomain.__init__(self,fileprefix=fileprefix,vertices=vertices,segments=segments,
                                               holes=holes,regions=regions,vertexFlags=vertexFlags,segmentFlags=segmentFlags,regionFlags=regionFlags,name=name,units=units)
        self.bathy=bathy
        self.bathyGridDim=bathyGridDim


class TriangulatedSurfaceDomain(D_base):
    """
    3D domains described by closed surfaces made up of triangular facets.
    """
    def __init__(self,vertices,triangles,
                 name="DefaultTriangulatedSurfaceDomain",units="m"):
        D_base.__init__(self,3,name,units)
        self.nodes=nodes
        self.triangles=triangles
    def readSTL(self,fileprefix):
        """
        Read the triangulated surface from a file in stereo lithography format.
        """
        raise UserWarning("STL input not implemented")
    def writePoly(self,fileprefix):
        """
        Write the triangulated surface in the poly format
        """
        raise UserWarning("Poly output not implemented")
    def writeAsymptote(self,fileprefix):
        """
        Write the triangulated surface in the Asymptote language.
        """
        raise UserWarning("Xdmf output not implemented")


class Mesh2DMDomain(D_base):
    """
    2D domains from ADH mesh files
    """
    def __init__(self,fileprefix):
        D_base.__init__(self,2,name=fileprefix)
        self.meshfile=fileprefix


class Mesh3DMDomain(D_base):
    """
    3D domains from ADH mesh files
    """
    def __init__(self,fileprefix):
        D_base.__init__(self,3,name=fileprefix)
        self.meshfile=fileprefix


class MeshHexDomain(D_base):
    """
    3D domains from Hex mesh files
    """
    def __init__(self,fileprefix):
        D_base.__init__(self,3,name=fileprefix)
        self.meshfile=fileprefix


class MeshTetgenDomain(D_base):
    """
    3D domains from tetgen mesh files
    """
    def __init__(self,fileprefix):
        D_base.__init__(self,3,name=fileprefix)
        self.meshfile=fileprefix


class PiecewiseLinearComplexDomain(D_base):
    """
    3D domains desribed by closed surfaces made up of general polygonal facets.
    """
    def __init__(self, fileprefix=None, vertices=None, facets=None,
                 facetHoles=None, holes=None, regions=None, vertexFlags=None, facetFlags=None,
                 regionFlags=None, regionConstraints=None, bc=None, name="PLCDomain", units="m"):
        """
        Read  the PLC domain from lists. If no flags are given, then default flags of 0 are assigned.
        """
        D_base.__init__(self, 3, name, units)
        if fileprefix:
            self.readPoly(fileprefix)
            self.regionLegend = {}
            self.regionConstraints = {}
            self.boundaryFlags = {}
        else:
            self.polyfile = None
            self.vertices = vertices or []
            self.facets = facets or []
            self.facetHoles = facetHoles or []
            self.holes = holes or []
            self.regions = regions or []
            self.vertexFlags = vertexFlags or []
            self.facetFlags = facetFlags or []
            self.regionFlags = regionFlags or []
            self.regionConstraints = regionConstraints or []
            self.regionLegend = {}
            self.boundaryFlags = {}
            import numpy as np
            self.update()

    def update(self):
        """
        Function to call when new vertices are added to the domain
        """
        if self.vertices:
            self.getBoundingBox()

    def readPoly(self,fileprefix):
        """
        Read in the PLC domain from a file in poly format.
        """
        self.polyfile = fileprefix
        self.name = fileprefix
        f = open(self.polyfile+".poly",'r')
        firstLine = f.readline().split()
        while len(firstLine) == 0 or firstLine[0][0] == '#':
            firstLine = f.readline().split()
        nVertices = int(firstLine[0])
        self.dim = int(firstLine[1])
        assert(self.dim == 3)
        nVertexAttributes = int(firstLine[2])
        self.vertexAttributes = [[] for j in range(nVertexAttributes)]
        hasVertexFlag = bool(int(firstLine[3]))
        self.vertices = []
        self.base = None
        self.vertexFlags = None
        if hasVertexFlag:
            self.vertexFlags = []
        for i in range(nVertices):
            line = f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            if not self.base:
                self.base = int(line[0])
            self.vertices.append([float(line[1]),float(line[2]),float(line[3])])
            for j in range(nVertexAttributes):
                self.vertexAttributes[j].append(float(line[4+j]))
            if hasVertexFlag:
                self.vertexFlags.append(int(line[4+nVertexAttributes]))
        facetLine = f.readline().split()
        while len(facetLine) == 0 or facetLine[0][0] == '#':
            facetLine = f.readline().split()
        nFacets = int(facetLine[0])
        hasFacetFlag = bool(int(facetLine[1]))
        self.facets = []
        self.facetHoles = []
        self.facetFlags = None
        if hasFacetFlag:
            self.facetFlags = []
        for i in range(nFacets):
            line = f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            self.facets.append([])
            self.facetHoles.append([])
            nPolygons = int(line[0])
            nHoles = int(line[1])
            if hasFacetFlag:
                self.facetFlags.append(int(line[2]))
            for j in range(nPolygons):
                line = f.readline().split()
                while len(line) == 0 or line[0][0] == '#':
                    line = f.readline().split()
                nSegments = int(line[0])
                self.facets[-1].append([int(line[1+k])-self.base for k in range(nSegments)])
            for j in range(nHoles):
                line = f.readline().split()
                while len(line) == 0 or line[0][0] == '#':
                    line = f.readline().split()
                self.facetHoles[-1].append([float(line[0]),
                                            float(line[1]),
                                            float(line[2])])
        holeLine = f.readline().split()
        while len(holeLine) == 0 or holeLine[0][0] == '#':
            holeLine = f.readline().split()
        nHoles = int(holeLine[0])
        self.holes = []
        for i in range(nHoles):
            line =  f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            self.holes.append([float(line[0]),float(line[1])])
        regionLine = f.readline().split()
        while len(regionLine) == 0 or regionLine[0][0] == '#':
            regionLine = f.readline().split()
        nRegions = int(regionLine[0])
        self.regions = []
        self.regionFlags = []
        self.areaConstraints = []
        for i in range(nRegions):
            line =  f.readline().split()
            while len(line) == 0 or line[0][0] == '#':
                line = f.readline().split()
            self.regions.append([float(line[1]),float(line[2]),float(line[3])])
            if len(line) > 4:
                self.regionFlags.append(int(line[4]))
            if len(line) > 5:
                self.areaConstraints.append(float(line[5]))
        self.getBoundingBox()
        f.close()

    def writePoly(self, fileprefix):
        """
        Write the PLC domain in the poly format.
        """
        if not self.polyfile:
            self.polyfile = fileprefix
            pf = open(fileprefix+'.poly','w')
            #write first line of poly file
            pf.write("%d %d %d %d \n" % (len(self.vertices), 3, 0, int(bool(self.vertexFlags))))
            #write the vertices
            for vN, v in enumerate(self.vertices):
                pf.write('%d %21.16e %21.16e %21.16e ' % (vN+1, v[0], v[1], v[2]))
                if self.vertexFlags :#write vertex flag if we have vertexFlags
                    pf.write('%d\n' % (self.vertexFlags[vN]))
                else:
                    pf.write('\n')
            #write the facets
            pf.write('%d %d \n' % (len(self.facets), int(bool(self.facetFlags))))
            for fN, f in enumerate(self.facets):
                if self.facetHoles:
                    nFacetHoles = len(self.facetHoles[fN])
                else:
                    nFacetHoles = 0
                if self.facetFlags:
                    pf.write('%d %d %d\n' % (len(f), nFacetHoles, self.facetFlags[fN]))
                else:
                    pf.write('%d %d\n' % (len(f), nFacetHoles))
                for segmentList in f:
                    pf.write(`len(segmentList)`+" ")
                    for vN in segmentList:
                        pf.write(`vN+1`+" ")
                    pf.write('\n')
                if self.facetHoles:
                    for hN, h in enumerate(self.facetHoles[fN]):
                        pf.write(`hN+1`+' %f %f %f\n' % h)
            if self.holes:
                pf.write('%d\n' % (len(self.holes),))
                for hN, h in enumerate(self.holes):
                    pf.write('%d %21.16e %21.16e %21.16e \n' % (hN+1,
                                                                h[0],
                                                                h[1],
                                                                h[2]))
            else:
                pf.write('%d\n' % (0,))
            if self.regions:
                pf.write('%d\n' % (len(self.regions),))
                for rN,r in enumerate(self.regions):
                    pf.write('%d %21.16e %21.16e %21.16e ' % (rN+1,
                                                              r[0],
                                                              r[1],
                                                              r[2]))
                    if self.regionFlags:
                        pf.write('%d ' % (self.regionFlags[rN]))
                    if self.regionConstraints:
                        pf.write('%21.16e' % (self.regionConstraints[rN]))
                    pf.write('\n')
            else:
                pf.write('%d\n' % (0,))
            pf.close()
        else:
            print "File already exists, not writing polyfile: " +`self.polyfile`

    def writePLY(self, fileprefix):
        """
        Write the PLC domain in the Stanford PLY format.
        """
        if True:#overwrite
            pf = open(fileprefix+'.ply','w')
            hasVertexFlags = int(bool(self.vertexFlags))
            hasFacetFlags=int(bool(self.facetFlags))
            hasFacetFlags=False
            hasVertexFlags=False # ?????
            pf.write("""ply
format ascii 1.0
comment author: Proteus
comment object: %s
""" % (self.name,))
            pf.write("""element vertex %d
property float x
property float y
property float z
""" % (len(self.vertices),))
            if hasVertexFlags:
                pf.write("property int flag\n")
            nSimpleFacets=0
            for f in self.facets:
                nSimpleFacets+=len(f)
            pf.write("element face %d \n" % (nSimpleFacets,))
            pf.write("property list uchar int vertex_index\n")
            if hasFacetFlags:
                pf.write("property int flag\n")
            pf.write("end_header\n")
            #write the vertices
            for vN,v in enumerate(self.vertices):
                pf.write('%21.16e %21.16e %21.16e' % (v[0],v[1],v[2]))
                if hasVertexFlags:
                    pf.write(" %d\n" % (self.vertexFlags[vN],))
                else:
                    pf.write("\n")
            #write the facets
            for fN,f in enumerate(self.facets):
                for segmentList in f:
                    pf.write(`len(segmentList)`)
                    for vN in segmentList:
                        pf.write(" "+`vN`)
                    if hasFacetFlags:
                        pf.write(" %d\n" % (self.facetFlags[fN],))
                    else:
                        pf.write("\n")
            pf.close()
        else:
            print "File already exists, not writing polyfile: " +`self.polyfile`

    def writeAsymptote(self, fileprefix):
        """
        Write a representation of the domain in the Asymptote vector graphics language.
        """
        f = open(fileprefix+'.asy', 'w')
        if not self.regionFlags:
            nM = 0
        else:
            nM = len(self.regionFlags)
#currentprojection=perspective(-2,-2,1,up=Z,target=(%(xc)f,%(yc)f,%(zc)f),showtarget=true,autoadjust=true,center=false);
        fileString = """import three;
import palette;
pen[] allPens=Wheel();
pen[] materialPens = new pen[%(nM)i];
for(int i=0;i< %(nM)i;++i)
  {
   int iPen = round(i*allPens.length/%(nM)i);
   materialPens[i] = allPens[iPen];
  }
currentprojection=FrontView;
unitsize(4.0 inches/%(Lx)f);
size(5 inches);
real Lx=%(Lx)f;
real Ly=%(Ly)f;
real Lz=%(Lz)f;
real offset=.0125Lx;
real x=%(x)f;
real y=%(y)f;
real z=%(z)f;
triple p0=(x,y,z);
triple p1=(x,y+Ly,z);
triple p2=(x+Lx,y+Ly,z);
triple p3=(x+Lx,y,z);
triple p4=(x,y,z+Lz);
triple p5=(x,y+Ly,z+Lz);
triple p6=(x+Lx,y+Ly,z+Lz);
triple p7=(x+Lx,y,z+Lz);
string strx="$%(Lx)2.2f\mbox{%(units)s}$";
string stry="$%(Ly)2.2f\mbox{%(units)s}$";
string strz="$%(Lz)2.2f\mbox{%(units)s}$";
draw(strx,(x,y-offset,z-offset)--(x+Lx,y-offset,z-offset),-(Y+Z),black,Bars3(Z),Arrows3);
draw(stry,(x-offset,y,z-offset)--(x-offset,y+Ly,z-offset),-(X+Z),black,Bars3(Z),Arrows3);
draw(strz,(x-offset,y-offset,z)--(x-offset,y-offset,z+Lz),-(X+Y),black,Bars3(X),Arrows3);
pen[] allPens=Wheel();
pen[] myPens = new pen[%(nFlags)i];
for(int i=0;i< %(nFlags)i;++i)
  {
   int iPen = round(i*allPens.length/%(nFlags)i);
   myPens[i] = allPens[iPen];
  }
""" % {'Lx':self.L[0],'Ly':self.L[1],'Lz':self.L[2],'x':self.x[0],'y':self.x[1],'z':self.x[2],'units':self.units,
       'xc':(self.x[0]+0.5*self.L[0]),
       'yc':(self.x[1]+0.5*self.L[1]),
       'zc':(self.x[2]+0.5*self.L[2]),
       'nM':nM,
       'nFlags':len(self.boundaryFlags)}
        f.write(fileString)
        for fN,facet in enumerate(self.facets):
            f.write(r"""surface s; path3 p;""")
            for vertexList in facet:
                if len(vertexList) > 1:
                    f.write("p = ((%f,%f,%f)" % tuple(self.vertices[vertexList[0]]))
                    for vN in vertexList[1:]:
                        f.write("--(%f,%f,%f)" % tuple(self.vertices[vN]))
                    f.write("--cycle);\n")
                    #f.write("draw(p);\n")
                    f.write("s.append(surface(p,planar=true));\n")
                #else:
                #    f.write("dot((%f,%f,%f));" % tuple(self.vertices[vertexList[0]]))
            #f.write("purge(10);\n")
        if self.regions and self.regionFlags :
            for region,regionFlag in zip(self.regions,self.regionFlags):
                f.write("dot((%f,%f,%f),materialPens[%i]);" % (region[0],region[1],region[2],regionFlag))
        f.close()
        f_bl = open(fileprefix+'_boundaryLegend.asy','w')
        fileString="""import palette;
pen[] allPens=Wheel();
unitsize(1.0 pt);
pen[] allPens=Wheel();
pen[] myPens = new pen[%(nFlags)i];
for(int i=0;i< %(nFlags)i;++i)
  {
   int iPen = round(i*allPens.length/%(nFlags)i);
   myPens[i] = allPens[iPen];
  }
dotfactor=12;
""" % {'nFlags':len(self.boundaryFlags)}
        f_bl.write(fileString)
        for bT,bN in self.boundaryFlags.iteritems():
            if bN != 0:
                f_bl.write("dot((12,%f),myPens[%i],L=Label(\"$%s$\",black));" % (bN*12,bN,bT))
        f_bl.close()

    def writeXdmf(self,fileprefix):
        pass




if __name__ == "__main__":
    import os
 #    r1d = RectangularDomain(L=[1.2])
#     r1d.writeAsymptote("r1d")
#     os.system("asy -V r1d")

#     r2d = RectangularDomain(L=[1.5,2.0],
#                             units='m')
#     r2d.writeAsymptote("r2d")
#     os.system("asy -V r2d")
#     r2d.writePoly("r2d")

#     r3d = RectangularDomain(L=[1.0,2.0,3.7],
#                             units='m')
#     r3d.writeAsymptote("r3d")
#     os.system("asy -V r3d")
#     r3d.writePoly("r3d")

#     plsg = PlanarStraightLineGraphDomain(vertices=[[0.0,0.0],
#                                                    [0.0,1.0],
#                                                    [1.0,0.0]],
#                                          segments=[[0,1],
#                                                    [1,2],
#                                                    [2,0]],
#                                          units='m')
#     plsg.writeAsymptote("pslg")
#     os.system("asy -V pslg")
#     plsg.writePoly("pslg")

#     plsgFromFile = PlanarStraightLineGraphDomain(units='m')
#     plsgFromFile.readPoly('pslg')
#     plsgFromFile.writeAsymptote('pslg2')
#     os.system("asy -V pslg2")

    plc = PiecewiseLinearComplexDomain(vertices=[[0.0,0.0,0.0],
                                                 [0.0,1.0,0.0],
                                                 [1.0,1.0,0.0],
                                                 [1.0,0.0,0.0],
                                                 [0.0,0.0,1.0],
                                                 [0.0,1.0,1.0],
                                                 [1.0,1.0,1.0],
                                                 [1.0,0.0,1.0],
                                                 [0.5,0.5,0.5],
                                                 [0.25,0.25,0.25],
                                                 [0.75,0.25,0.25]],
                                       facets=[[[0,1,2,3]],
                                               [[0,4,7,3]],
                                               [[0,4,5,1]],
                                               [[4,5,6,7]],
                                               [[3,7,6,2]],
                                               [[2,6,5,1]],
                                               [[8]],
                                               [[9,10]]],
                                       units='m')
    plc.writeAsymptote("plc")
    os.system("asy -V plc")
    plc.writePoly("plc")
    plc.writeGeo("geo")

#    plcFromFile = PiecewiseLinearComplexDomain(units='m')
#    plcFromFile.readPoly('cylinder3d')
#    plcFromFile.writeAsymptote('cylinder3d')
#    os.system("asy -V cylinder3d")
#    plcFromFile = PiecewiseLinearComplexDomain(units='m')
#    plcFromFile.readPoly('Container06.GF')
#    plcFromFile.writeAsymptote('Container06.GF')
#    os.system("asy -V Container06.GF")
