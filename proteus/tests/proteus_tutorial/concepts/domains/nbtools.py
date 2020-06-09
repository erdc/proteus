from pythreejs import *
from IPython.display import set_matplotlib_formats,display
from matplotlib import pyplot, collections
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
from matplotlib import pyplot as plt
import pylab
from proteus import SpatialTools as st
import numpy as np

def plot_domain(domain, elev=0, azim=0, flags=None):
    bc = domain.shape_list[0].BC_class
    st._assembleGeometry(domain, bc)
    if domain.nd == 2:
        lines = []
        cmap = pylab.get_cmap("hsv")
        c = []
        shape_nb = float(len(domain.shape_list))
        annotate_coords = []
        annotate_label = []
        start_global_flag = 0
        if flags=='global': shapes = [domain];
        if flags=='local' or flags is None: shapes = domain.shape_list
        for i, shape in enumerate(shapes):
            for s,sF in zip(shape.segments,shape.segmentFlags):
                lines.append([shape.vertices[s[0]], shape.vertices[s[1]]])
                c.append(cmap(float(i/shape_nb)))
                annotate_coords += [np.array((np.array(shape.vertices[s[0]])+np.array(shape.vertices[s[1]]))/2.).tolist()]
                annotate_label += [sF]
        lc = collections.LineCollection(lines,colors=c,linewidths=3)
        fig, ax = pyplot.subplots()
        if flags is not None:
            for i, label in enumerate(annotate_label):
                ax.annotate(str(label), xy=annotate_coords[i], xycoords='data',
                #bbox=dict(boxstyle="round", fc="0.8"))
                )
        ax.add_collection(lc)
        ax.margins(0.1)
        ax.set_aspect('equal')
        #pylab.savefig("pslg.pdf")
    if domain.nd == 3:
        ax = a3.Axes3D(pl.figure(),elev=elev,azim=azim)
        cmap = pylab.get_cmap("hsv")
        fN_max = float(max(domain.facetFlags))
        verts=[]
        c=[]
        for f,fN in zip(domain.facets,domain.facetFlags):
            verts.append([domain.vertices[vN] for vN in f[0]])
            #c.append(cmap(fN/fN_max))
        ply = a3.art3d.Poly3DCollection(verts)
        #ply.set_facecolors(c)
        ax.add_collection3d(ply)
        ax.margins(0.1,0.1,0.1)
        #ax.set_aspect('equal')
        ax.set_xlim(domain.x[0],max(domain.L))
        ax.set_ylim(domain.x[1],max(domain.L))
        ax.set_zlim(domain.x[2],max(domain.L))
        pylab.savefig("plc.pdf")

def plot_js(domain):
    bc = domain.shape_list[0].BC_class
    st._assembleGeometry(domain, bc)
    faces3=[]
    facesn=[]
    vertices=[]
    va = np.array(domain.vertices)
    cg = va.mean(0)
    for v in domain.vertices:
        vertices+=(np.array(v)-cg).tolist()
    for facet in domain.facets:
        for face in facet:
            facesn.append(face)
            for i in range(len(face)-2):
                faces3.append([face[0],face[i+1],face[i+2]])
    p  = FaceGeometry(face3=[],face4=[],facen=facesn,vertices=vertices)
    mat = BasicMaterial(side='DoubleSide',color='red')
    m = Mesh(geometry=p, material=mat)
    scene = Scene(children=[m, AmbientLight(color=0x777777)])
    scene.background=0x777777
    c = PerspectiveCamera(position=[0,domain.x[1]+1.5*domain.L[1],domain.x[2]+1.5*domain.L[2]], up=[0,0,1], 
                          children=[DirectionalLight(color=0x777777, position=[3,5,1], intensity=0.6)])
    renderer = Renderer(camera=c, scene = scene, controls=[OrbitControls(controlling=c)],background="grey")
    display(renderer)