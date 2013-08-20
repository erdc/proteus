from proteus.iproteus import *
from proteus import Comm
from proteus.Domain import InterpolatedBathymetryDomain
from proteus.MeshTools import InterpolatedBathymetryMesh
from proteus.Archiver import XdmfArchive
import xml.etree.ElementTree as ElementTree

comm = Comm.init()
Profiling.logLevel=7
Profiling.verbose=True

def setupStepGauss():
    import numpy as np
    from math import sin,cos,pi,sqrt,exp
    #set up a fake LiDAR point set
    nPoints_x = nPoints_y = 51
    delta_x = 2.0/float(nPoints_x-1)
    delta_y = 2.0/float(nPoints_y-1)
    bathy = np.zeros((nPoints_x*nPoints_y,3),'d')
    for i in range(nPoints_y):
        for j in range(nPoints_x):
            x = -0.5+j*delta_x
            y = -0.5+i*delta_y
            #z = 1.0
            if y > x:
                z = 1.0
            else:
                if y < x - 0.25:
                    r = sqrt((y - 0.25)**2 + (x - 0.8)**2)
                    z = exp(-50.0*r**2)
                else:
                    z = 0.0
            #z = y
            #z += sin(2.0*pi*x)*cos(2.0*pi*y)
            bathy[i*nPoints_x+j,0] = x
            bathy[i*nPoints_x+j,1] = y
            bathy[i*nPoints_x+j,2] = z
    domain = InterpolatedBathymetryDomain(vertices=[[0.0,0.0],
                                                    [0.0,1.0],
                                                    [0.5,1.5],
                                                    [1.0,1.0],
                                                    [1.5,-0.5]],
                                          vertexFlags=[1,2,3,2,1],
                                          segments=[[0,1],
                                                    [1,2],
                                                    [2,3],
                                                    [3,4],
                                                    [4,0]],
                                          segmentFlags=[1,2,3,3,1],
                                          regions=[(0.5,0.5)],
                                          regionFlags=[1],
                                          name="interpolatedBathySimpleTest",
                                          units='m',
                                          tol=max(delta_x,delta_y),
                                          bathy = bathy)
    domain.writePoly(domain.name)
    return domain

def test_L1():
    domain = setupStepGauss()
    mesh = InterpolatedBathymetryMesh(domain,triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),maxLevels=25,maxNodes=50000,normType="L1")
    archive = XdmfArchive(dataDir='.',filename="interpolatedBathySimpleTest_L1_")
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
    archive.close()

def test_L2():
    domain = setupStepGauss()
    mesh = InterpolatedBathymetryMesh(domain,triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),maxLevels=25,maxNodes=50000,normType="L2")
    archive = XdmfArchive(dataDir='.',filename="interpolatedBathySimpleTest_L2_")
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
    archive.close()

def test_Linfty():
    domain = setupStepGauss()
    mesh = InterpolatedBathymetryMesh(domain,triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),maxLevels=25,maxNodes=50000,normType="Linfty")
    archive = XdmfArchive(dataDir='.',filename="interpolatedBathySimpleTest_Linfty_")
    archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
    mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
    archive.close()

if __name__ == '__main__':
    test_L1()
    test_L2()
    test_Linfty()
