import os

from proteus.Domain import InterpolatedBathymetryDomain
from proteus.MeshTools import InterpolatedBathymetryMesh
from proteus.Archiver import XdmfArchive

import xml.etree.ElementTree as ElementTree
import pytest

@pytest.mark.MeshTools
@pytest.mark.Archiver
@pytest.mark.Domain
class TestInterpolatedBathy():
    """ Runs a set of tests for Interpolated Bathymetry"""
    
    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        self.aux_names = []

    def teardown_method(self,method):
        filenames = []
        for aux_name in self.aux_names:
            filenames.extend([aux_name+'.'+post for post in ['xmf','h5','2dm']])
            filenames.extend([aux_name+'0.'+post for post in ['xmf','h5','2dm']])
        filenames.extend(['tetgen'+'.'+post for post in ['ele','node','face']])
        filenames.extend(['proteus_default.log','interpolatedBathySimpleTest.poly'])
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError, e:
                    print ("Error: %s - %s." %(e.filename,e.strerror))
            else:
                pass

    def setupStepGauss(self):
        import numpy as np
        from math import sin,cos,pi,sqrt,exp
        #set up a fake LiDAR point set
        nPoints_x = nPoints_y = 21
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
                                              bathy = bathy,
                                              bathyGridDim = (nPoints_y,nPoints_x))
        domain.writePoly(domain.name)
        return domain

    def test_L1(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="points",
                                          bathyAssignmentScheme="localAveraging",errorNormType="L1")
        outfile = "interpolatedBathySimpleTest_L1_"
        archive = XdmfArchive(dataDir='.',filename = outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_L2(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="points",
                                          bathyAssignmentScheme="localAveraging",
                                          errorNormType="L2")
        outfile = "interpolatedBathySimpleTest_L2_"
        archive = XdmfArchive(dataDir='.',filename = outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_Linfty(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="points",
                                          bathyAssignmentScheme="localAveraging",
                                          errorNormType="Linfty")
        outfile = "interpolatedBathySimpleTest_Linfty_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)
        
    def test_L1_interp(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="points",
                                          bathyAssignmentScheme="interpolation",
                                          errorNormType="L1")
        outfile = "interpolatedBathySimpleTest_L1_interp_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)
        
    def test_L2_interp(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="points",
                                          bathyAssignmentScheme="interpolation",
                                          errorNormType="L2")
        outfile = "interpolatedBathySimpleTest_L2_interp_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_Linfty_interp(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="points",
                                          bathyAssignmentScheme="interpolation",
                                          errorNormType="Linfty")
        outfile = "interpolatedBathySimpleTest_Linfty_interp_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_L1_grid(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="grid",
                                          bathyAssignmentScheme="localAveraging",
                                          errorNormType="L1")
        outfile = "interpolatedBathySimpleTest_grid_L1_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_L2_grid(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="grid",
                                          bathyAssignmentScheme="localAveraging",
                                          errorNormType="L2")
        outfile = "interpolatedBathySimpleTest_grid_L2_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_Linfty_grid(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="grid",
                                          bathyAssignmentScheme="localAveraging",
                                          errorNormType="Linfty")
        outfile = "interpolatedBathySimpleTest_grid_Linfty_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_L1_interp_grid(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="grid",
                                          bathyAssignmentScheme="interpolation",
                                          errorNormType="L1")
        outfile = "interpolatedBathySimpleTest_grid_L1_interp_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_L2_interp_grid(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="grid",
                                          bathyAssignmentScheme="interpolation",
                                          errorNormType="L2")
        outfile = "interpolatedBathySimpleTest_grid_L2_interp_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        self.aux_names.append(outfile)

    def test_Linfty_interp_grid(self):
        domain = self.setupStepGauss()
        mesh = InterpolatedBathymetryMesh(domain,
                                          triangleOptions="gVApq30Dena%8.8f" % (0.5**3,),
                                          atol=1.0e-1,
                                          rtol=1.0e-1,
                                          maxLevels=25,
                                          maxNodes=50000,
                                          bathyType="grid",
                                          bathyAssignmentScheme="interpolation",
                                          errorNormType="Linfty")
        outfile = "interpolatedBathySimpleTest_grid_Linfty_interp_"
        archive = XdmfArchive(dataDir='.',filename=outfile, global_sync=False)
        archive.domain = ElementTree.SubElement(archive.tree.getroot(),"Domain")
        mesh.meshList[-1].writeMeshXdmf(ar=archive,init=True)
        archive.sync(); archive.close()
        mesh.meshList[-1].writeMeshADH("interpolatedBathySimpleTest_grid_Linfty_interp_")
        self.aux_names.append(outfile)


if __name__ == '__main__':
    pass
