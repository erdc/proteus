from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = (5.0*4.0/4.0, 4.0)
import numpy as np
from IPython.display import set_matplotlib_formats,display
set_matplotlib_formats('svg')
def plot_rectangular_domain(domain):
    from matplotlib import collections
    cmap = plt.get_cmap("hsv")
    lines = [[domain.x,(domain.x[0],domain.x[1]+domain.L[1])], 
             [(domain.x[0],domain.x[1]+domain.L[1]),(domain.x[0]+domain.L[0],domain.x[1]+domain.L[1])],
             [(domain.x[0]+domain.L[0],domain.x[1]+domain.L[1]),(domain.x[0]+domain.L[0],domain.x[1])],
             [(domain.x[0]+domain.L[0],domain.x[1]),domain.x]]
    lc = collections.LineCollection(lines,linewidths=3)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.margins(0.1)
    ax.set_aspect('equal')
    plt.savefig("rectangle.pdf")
def plot_pslg_domain(domain):
    from matplotlib import collections
    lines = []
    cmap = plt.get_cmap("hsv")
    c = []
    text_points = []
    sF_max = float(max(domain.segmentFlags))
    for s,sF in zip(domain.segments,domain.segmentFlags):
        lines.append([domain.vertices[s[0]],domain.vertices[s[1]]])
        text_points.append(((domain.vertices[s[0]][0] + domain.vertices[s[1]][0])/2.,
                            (domain.vertices[s[0]][1] + domain.vertices[s[1]][1])/2.))
        c.append(cmap(float(sF/sF_max)))
    lc = collections.LineCollection(lines,colors=c,linewidths=3)
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    for p,sF in zip(text_points,domain.segmentFlags):
        ax.text(p[0],p[1],str(sF))
    ax.margins(0.1)
    ax.set_aspect('equal')
    plt.savefig("pslg.pdf")
def plot_plc_domain(domain):
    import mpl_toolkits.mplot3d as a3
    import matplotlib.colors as colors
    import matplotlib.patches as mpatches
    ax = a3.Axes3D(plt.figure())
    cmap = plt.get_cmap("hsv")
    fN_max = float(max(domain.facetFlags))
    verts=[]
    c=[]
    flags={}
    for f,fN in zip(domain.facets,domain.facetFlags):
        verts.append([domain.vertices[vN] for vN in f[0]])
        c.append(cmap(fN/fN_max))
        flags[fN] = mpatches.Patch(color=c[-1], label=str(fN))
    ply = a3.art3d.Poly3DCollection(verts)
    ply.set_facecolors(c)
    ax.add_collection3d(ply)
    ax.margins(0.1,0.1,0.1)
    #ax.set_aspect('equal')
    ax.set_xlim(domain.x[0],domain.L[0])
    ax.set_ylim(domain.x[1],domain.L[1])
    ax.set_zlim(domain.x[2],domain.L[2])
    plt.legend(handles=flags.values())
    plt.savefig("plc.pdf")