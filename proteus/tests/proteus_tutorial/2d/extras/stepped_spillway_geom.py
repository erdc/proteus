
# coding: utf-8

# In[218]:


import numpy as np
#import matplotlib.pyplot as plt
from proteus import SpatialTools as st
import os
from proteus import Domain

def get_yz(filename):
    data=st.getInfoFromSTL(str(filename))
    x=data[0][0][0]
    yz=[]
    for i in range(len(data[0])):
        if data[0][i][0] >= 0.9999*x and data[0][i][0] <= 1.0001*x:
            yz.append([0.0,data[0][i][1],data[0][i][2]])
    return yz

def yz2xy(yz):
    xy=[]
    for i in range(len(yz)):
        xy.append([yz[i][1],yz[i][2]])
    xy.sort(key=lambda tup: tup[0])
    return xy

def clipds(vrts):
    endl=vrts[-3][0]
    newl=[]
    for i in range(len(vrts)):
        if (vrts[i][0]> 426) and (vrts[i][0]<endl):
            pass
        else:
            newl.append(vrts[i])
    return newl

def clipus(vrts):
    begl=vrts[0][0]
    newl=[]
    for i in range(len(vrts)):
        if (vrts[i][0]> begl) and (vrts[i][0]<303):
            pass
        else:
            newl.append(vrts[i])
    return newl

#def plots(fig):
#    plt.figure(figsize=(40,10))
#    plt.scatter(np.array(fig)[:,0],np.array(fig)[:,1])
#    plt.show()

def v_flags(vertices,boundaryTags):
    vertexFlags=[]
    for i in range(len(vertices)):
        if i < (len(vertices)-3):
                vertexFlags.append(boundaryTags['bottom'])
        elif i == (len(vertices)-3):
                vertexFlags.append(boundaryTags['outflow'])
        elif i == (len(vertices)-2):
                vertexFlags.append(boundaryTags['top'])
        elif i == (len(vertices)-1):
                vertexFlags.append(boundaryTags['inflow'])
    return vertexFlags

def segs(vertices):
    segments=[]
    for i in range(len(vertices)):
        a=[]
        if i < len(vertices)-1:
            segments.append([i,i+1])
        else:
            segments.append([i,0])
    return segments

def s_flags(segments,boundaryTags):
    segmentFlags=[]
    for i in range(len(segments)):
        if i < (len(segments)-3):
                segmentFlags.append(boundaryTags['bottom'])
        elif i == (len(segments)-3):
                segmentFlags.append(boundaryTags['outflow'])
        elif i == (len(segments)-2):
                segmentFlags.append(boundaryTags['top'])
        elif i == (len(segments)-1):
                segmentFlags.append(boundaryTags['inflow'])
    return segmentFlags

def geom_transform(vertices):
    a=min(np.array(vertices)[:,0])
    b=min(np.array(vertices)[:,1])
    for i in range(len(vertices)):
        vertices[i]=[vertices[i][0]-a,vertices[i][1]-b]
    return vertices

def top(vertices,top):
    vertices[-1][1]=top
    vertices[-2][1]=top
    return vertices

def ups_len(vertices,upslen):
    if upslen==None:
        return vertices
    else:
        a=vertices[1][0]-vertices[0][0]
        for i in range(len(vertices)):
            if (i > 0 and i < len(vertices)-1):
                vertices[i][0]=(vertices[i][0]-a+upslen)
            else:
                pass
        return vertices

def dwns_len(vertices,dwnslen):
    if dwnslen ==None:
        return vertices
    else:
        a=vertices[-3][0]-vertices[-4][0]
        for i in range(len(vertices)):
            if (i < len(vertices)-1 and i > len(vertices)-4):
                vertices[i][0]=(vertices[i][0]-a+dwnslen)
            else:
                pass
        return vertices

def deldup(vertices):
    a=[vertices[0]]
    for i in range(len(vertices)):
        if i == 0 or vertices[i]==vertices[i-1]:
            pass
        else:
            a.append(vertices[i])
    return a

def makecsv():
    a=get_yz('bed1.stl')
    b=get_yz('bed2.stl')
    c=get_yz('conc.stl')
    d=a+b+c

    np.savetxt("full_bed.csv",d, delimiter=",")

    f=[]
    e=yz2xy(d)
    f=e
    f.append([max(np.array(e)[:,0]),max(np.array(get_yz('outlet.stl'))[:,2])])
    f.append([min(np.array(f)[:,0]),max(np.array(get_yz('outlet.stl'))[:,2])])
    g=f
    h=clipus(clipds(g))
    i=deldup(h)
    np.savetxt("domain.csv",f, delimiter=",")
#    np.savetxt("domain_clip.csv",h, delimiter=",")
    np.savetxt("domain_clip.csv",i, delimiter=",")
    boundaries=['bottom','outflow','top','inflow']
    boundaryTags=dict([(key,i+1) for (i,key) in enumerate(boundaries)])

    vertices=np.genfromtxt("domain_clip.csv", delimiter=",").tolist()
    vertexFlags=v_flags(vertices,boundaryTags)
    segments=segs(vertices)
    segmentFlags=s_flags(segments,boundaryTags)

