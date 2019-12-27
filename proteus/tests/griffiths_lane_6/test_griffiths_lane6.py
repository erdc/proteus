#!/usr/bin/env python
"""
Test module for linear boundary value problems (serial)

This module solves equations of the form

.. math::

  \nabla \cdot \left( a(x) \nabla u \right) = f(x)

"""
from __future__ import print_function
from builtins import object
from proteus.iproteus import *
import os
from past.utils import old_div
import numpy as np
import tables
from . import re_gl_6_3d_p
from . import re_gl_6_3d_n
from . import sm_gl_6_3d_p
from . import sm_gl_6_3d_n

#BASELINE L2 norm value used for test comparison, expected.root.pressure_head_t1
#richards_L2_norm_pressure_baseline = 1192.0435669518072
richards_L2_norm_pressure_baseline = 1113.0020868267397

#elastoplastic expected.root.displacement_t1
#elastoplastic_L2_norm_displacement_baseline_x = 3.96325533e-01 
#elastoplastic_L2_norm_displacement_baseline_y = 1.33090951e-05
#elastoplastic_L2_norm_displacement_baseline_z = 1.85643753e+00
elastoplastic_L2_norm_displacement_baseline_x = 1.56320197e-01
elastoplastic_L2_norm_displacement_baseline_y = 1.66025679e-07
elastoplastic_L2_norm_displacement_baseline_z = 1.52434107e+00


from proteus import Quadrature
from proteus import MeshTools
from proteus import FemTools

def getDummyMesh(IEN,nodeIDs):
    dummyMesh = MeshTools.Mesh()
    dummyMesh.nodeOffsets_subdomain_owned = [0,0]
    dummyMesh.globalMesh = dummyMesh
    dummyMesh.max_nNodeNeighbors_node = 0
    dummyMesh.nElements_global = IEN.shape[0]
    dummyMesh.dim = 3

    dummyMesh.nodeArray=nodeIDs
    dummyMesh.elementNodesArray=IEN

    return dummyMesh


def get_L2_norm(h5file,field):

    IEN = np.array(h5file.root.elementsSpatial_Domain1)
    nodeIDs = np.array(h5file.root.Mesh_Spatial_Domain_1)
    nodeCoords = np.array(h5file.root.nodesSpatial_Domain1)

    dummyMesh = getDummyMesh(IEN,nodeIDs)
    
    elementQuadrature = Quadrature.SimplexGaussQuadrature(dummyMesh.dim, 2)
    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh=dummyMesh,nd=3)

    L2_norm_total = 0

    for eID,ele in enumerate(IEN):
        vertices = nodeCoords[ele]
        
        volumeMatrix = np.vstack([vertices.transpose(),[1,1,1,1]])
        volume = abs(np.linalg.det(volumeMatrix)/6.0)
        localBasis= femSpace.getBasisValuesRef(np.array(elementQuadrature.points))
        #jacobian
        J = volume/(old_div(1.0,6.0))
        scalar = np.dot(field[dummyMesh.elementNodesArray[eID]],localBasis)
        L2_norm = np.dot(np.square(scalar),np.array(elementQuadrature.weights))
        L2_norm = L2_norm*J
        L2_norm_total+=L2_norm

    L2_norm_total = np.sqrt(L2_norm_total)

    return L2_norm_total

def get_L2_vectorNorm(h5file,field):

    IEN = np.array(h5file.root.elementsSpatial_Domain1)
    nodeIDs = np.array(h5file.root.Mesh_Spatial_Domain_1)
    nodeCoords = np.array(h5file.root.nodesSpatial_Domain1)
    dummyMesh = getDummyMesh(IEN,nodeIDs)
    
    elementQuadrature = Quadrature.SimplexGaussQuadrature(dummyMesh.dim, 2)
    femSpace = FemTools.C0_AffineLinearOnSimplexWithNodalBasis(mesh=dummyMesh,nd=3)

    L2_norm_total = np.array([0.0,0.0,0.0])

    for eID,ele in enumerate(IEN):
        vertices = nodeCoords[ele]
        
        volumeMatrix = np.vstack([vertices.transpose(),[1,1,1,1]])
        volume = abs(np.linalg.det(volumeMatrix)/6.0)
        localBasis=femSpace.getBasisValuesRef(np.array(elementQuadrature.points))
        J = volume/(old_div(1.0,6.0))
        for i in range(dummyMesh.dim):
            vectorComp = np.dot(field[dummyMesh.elementNodesArray[eID],i],localBasis)
            L2_norm = np.dot(np.square(vectorComp),np.array(elementQuadrature.weights))
            L2_norm = L2_norm*J
            L2_norm_total[i]+=L2_norm

    L2_norm_total = np.sqrt(L2_norm_total)
    return L2_norm_total


class TestRichards(object):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        self.aux_names = []
        self.meshdir = os.path.dirname(os.path.abspath(__file__))
        self._scriptdir = os.path.dirname(os.path.abspath(__file__))
        
    def teardown_method(self,method):
        filenames = []
        #for aux_name in self.aux_names:
        #    filenames.extend([aux_name+'.'+ext for ext in ['h5','xmf']])
        filenames.append('proteus.log')
        for f in filenames:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError as e:
                    print ("Error: %s - %s" %(e.filename,e.strerror))
            else:
                pass
            
    def test_richards(self,use_strong_constraints=False):
        pList = [re_gl_6_3d_p]
        nList = [re_gl_6_3d_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "richards"
        reload(default_s)
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution(so.name)
        self.aux_names.append(so.name)
        # COMPARE VS SAVED FILES #
        actual = tables.open_file(so.name+'.h5','r')
        L2_norm = get_L2_norm(actual,actual.root.pressure_head_t1)
        np.testing.assert_almost_equal(L2_norm,richards_L2_norm_pressure_baseline)
        actual.close()
        del ns
        
    def test_elastoplastic(self,use_strong_constraints=False):
        pList = [sm_gl_6_3d_p]
        nList = [sm_gl_6_3d_n]
        reload(default_so)
        so = default_so
        so.name = pList[0].name = "elastoplastic"
        reload(default_s)
        so.sList=[default_s]
        opts.logLevel=7
        opts.verbose=True
        opts.profile=True
        opts.gatherArchive=True
        try:
            ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        except:
            assert 0, "Failed at NS_base"
        ns.calculateSolution(so.name)
        self.aux_names.append(so.name)
        # COMPARE VS SAVED FILES #
        expected_path = so.name+'_expected.h5'
        expected = tables.open_file(os.path.join(self._scriptdir,expected_path))
        actual = tables.open_file(so.name+'.h5','r')
        L2_norm = get_L2_vectorNorm(actual,actual.root.displacement_t1)
        baseline = np.array([elastoplastic_L2_norm_displacement_baseline_x, elastoplastic_L2_norm_displacement_baseline_y, elastoplastic_L2_norm_displacement_baseline_z])
        np.testing.assert_almost_equal(L2_norm,baseline)
        actual.close()
        del ns
        
if __name__ == '__main__':
    pass
