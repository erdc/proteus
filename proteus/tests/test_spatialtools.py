"""
Testing module for Domain.py, Shape.py, BC.py
Work in progress
"""
import unittest
import numpy.testing as npt
import numpy as np
from nose.tools import eq_
from proteus import Domain
from proteus import SpatialTools as st
from proteus.mprans import SpatialTools as strans


def create_domain2D():
    return Domain.PlanarStraightLineGraphDomain()

def create_domain3D():
    return Domain.PiecewiseLinearComplexDomain()

def create_rectangle(domain, dim=(0., 0.), coords=(0., 0.), folder=None):
    if folder is None:
        return st.Rectangle(domain, dim, coords)
    elif folder == 'mprans':
        return strans.Rectangle(domain, dim, coords)

def create_cuboid(domain, dim=(0., 0., 0.), coords=(0., 0., 0.), folder=None):
    if folder is None:
        return st.Cuboid(domain, dim, coords)
    elif folder == 'mprans':
        return strans.Cuboid(domain, dim, coords)

def create_custom2D(domain, folder=None):
    domain2D = domain
    bt2D = {'bottom': 1, 'right': 2, 'left': 3, 'top': 4}
    vertices2D = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]]
    vertexFlags2D = [bt2D['bottom'], bt2D['bottom'], bt2D['top'], bt2D['top']]
    segments2D = [[0, 1], [1, 2], [2, 3], [3, 0]]
    segmentFlags2D = [bt2D['bottom'], bt2D['right'], bt2D['top'], bt2D['left']]
    if folder is None:
        mod = st
    elif folder == 'mprans':
        mod = strans
    custom2D = mod.CustomShape(domain=domain2D, vertices=vertices2D,
                               vertexFlags=vertexFlags2D, segments=segments2D,
                               segmentFlags=segmentFlags2D, boundaryTags=bt2D)
    return custom2D

def create_custom3D(domain, folder=None):
    domain3D = domain
    bt3D = {'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5,
            'top':6}
    vertices3D = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                  [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]]
    vertexFlags3D = [bt3D['left'], bt3D['right'], bt3D['right'], bt3D['left'],
                     bt3D['left'], bt3D['right'], bt3D['right'], bt3D['left']]
    facets3D = [[[0, 1, 2, 3]], [[0, 1, 5, 4]], [[1, 2, 6, 5]], [[2, 3, 7, 6]],
                [[3, 0, 4, 7]], [[4, 5, 6, 7]]]
    facetFlags3D = [bt3D['bottom'], bt3D['front'], bt3D['right'], bt3D['back'],
                    bt3D['left'], bt3D['top']]
    if folder is None:
        mod = st
    elif folder == 'mprans':
        mod = strans
    custom3D = mod.CustomShape(domain=domain3D, vertices=vertices3D,
                               vertexFlags=vertexFlags3D, facets=facets3D,
                               facetFlags=facetFlags3D, boundaryTags=bt3D)
    return custom3D

def create_tank2D(domain, dim=(0., 0.), coords=(0., 0.)):
    return strans.Tank2D(domain, dim, coords)

def create_tank3D(domain, dim=(0., 0., 0.), coords=(0., 0., 0.)):
    return strans.Tank3D(domain, dim, coords)


class TestShapeDomainBuilding(unittest.TestCase):

    def test_create_shapes(self):
        """
        Testing if shapes can be created
        """
        domain2D = create_domain2D()
        domain3D = create_domain3D()
        rectangle = create_rectangle(domain2D)
        rectangleRANS = create_rectangle(domain2D, folder='mprans')
        cuboid = st.Cuboid(domain3D)
        cuboidRANS = create_cuboid(domain3D, folder='mprans')
        tand2D = create_tank2D(domain2D)
        tank3D = create_tank3D(domain3D)
        custom2D = create_custom2D(domain2D) 
        custom2DRANS = create_custom2D(domain2D, folder='mprans') 
        custom3D = create_custom3D(domain3D) 
        custom3DRANS = create_custom3D(domain3D, folder='mprans') 

    def test_assemble_domain(self):
        """
        Testing the assembleDomain() for different domains with multiple shapes
        """
        nb_shapes = 10
        domain2D = create_domain2D()
        domain2DRANS = create_domain2D()
        domain3D = create_domain3D()
        domain3DRANS = create_domain3D()
        dim2D = [1., 1.]
        dim3D = [1., 1., 1.]
        coords2D = np.array([0.5, 0.5])
        coords3D = np.array([0.5, 0.5, 0.5])
        nb_bc2D = 0
        nb_bc2DRANS = 0
        nb_bc3D = 0
        nb_bc3DRANS = 0
        for shape in range(nb_shapes):
            coords2D += 1.5
            a = create_rectangle(domain2D, dim=dim2D, coords=coords2D)
            nb_bc2D += len(a.BC_list)
            a = create_cuboid(domain3D, dim=dim3D, coords=coords3D)
            nb_bc3D += len(a.BC_list)
            a = create_rectangle(domain2DRANS, dim=dim2D, coords=coords2D,
                                folder='mprans')
            nb_bc2DRANS += len(a.BC_list)
            a = create_cuboid(domain3DRANS, dim=dim3D, coords=coords3D,
                            folder='mprans')
            nb_bc3DRANS += len(a.BC_list)
        a = create_tank2D(domain2DRANS, dim=[50., 50.], coords=[25., 25.])
        nb_bc2DRANS += len(a.BC_list)
        a = create_tank3D(domain3DRANS, dim=[50., 50., 50.],
                          coords=[25., 25., 25.])
        nb_bc3DRANS += len(a.BC_list)
        st.assembleDomain(domain2D)
        st.assembleDomain(domain3D)
        strans.assembleDomain(domain2DRANS)
        strans.assembleDomain(domain3DRANS)
        
        # check that each domain has the right number of shapes
        npt.assert_equal(len(domain2D.shape_list), nb_shapes)
        npt.assert_equal(len(domain3D.shape_list), nb_shapes)
        npt.assert_equal(len(domain2DRANS.shape_list), nb_shapes+1)
        npt.assert_equal(len(domain3DRANS.shape_list), nb_shapes+1)

        # check that the number of boundary conditions is right
        npt.assert_equal(len(domain2D.bc), nb_bc2D+1)
        npt.assert_equal(len(domain3D.bc), nb_bc3D+1)
        npt.assert_equal(len(domain2DRANS.bc), nb_bc2DRANS+1)
        npt.assert_equal(len(domain3DRANS.bc), nb_bc3DRANS+1)

class TestFlags(unittest.TestCase):

    def test_BC_flags(self):
        """
        Testing the flags of shapes and their in their domain
        """
        nb_shapes = 3
        domain2D = create_domain2D()
        domain3D = create_domain3D()
        domain2DRANS = create_domain2D()
        domain3DRANS = create_domain3D()
        flags_v2D = []
        flags_s2D = []
        flags_v3D = []
        flags_f3D = []
        flags_v2DRANS = []
        flags_s2DRANS = []
        flags_v3DRANS = []
        flags_f3DRANS = []
        maxf = 0
        for i in range(nb_shapes):
            # 2D
            a = create_custom2D(domain2D)
            if flags_v2D:
                maxf = np.max([np.max(flags_v2D), np.max(flags_s2D)])
            flags_v2D += (a.vertexFlags+maxf).tolist()
            flags_s2D += (a.segmentFlags+maxf).tolist()
            # 3D
            a = create_custom3D(domain3D)
            if flags_v3D:
                maxf = np.max([np.max(flags_v3D), np.max(flags_f3D)])
            flags_v3D += (a.vertexFlags+maxf).tolist()
            flags_f3D += (a.facetFlags+maxf).tolist()
            # 2D RANS
            a = create_custom2D(domain2DRANS, folder='mprans')
            if flags_v2DRANS:
                maxf = np.max([np.max(flags_v2DRANS), np.max(flags_s2DRANS)])
            flags_v2DRANS += (a.vertexFlags+maxf).tolist()
            flags_s2DRANS += (a.segmentFlags+maxf).tolist()
            # 3D RANS
            a = create_custom3D(domain3DRANS, folder='mprans')
            if flags_v3DRANS:
                maxf = np.max([np.max(flags_v3DRANS), np.max(flags_f3DRANS)])
            flags_v3DRANS += (a.vertexFlags+maxf).tolist()
            flags_f3DRANS += (a.facetFlags+maxf).tolist()
        st.assembleDomain(domain2D)
        st.assembleDomain(domain3D)
        strans.assembleDomain(domain2DRANS)
        strans.assembleDomain(domain3DRANS)

        # test flags
        npt.assert_equal(domain2D.vertexFlags, flags_v2D)
        npt.assert_equal(domain2D.segmentFlags, flags_s2D)
        npt.assert_equal(domain3D.vertexFlags, flags_v3D)
        npt.assert_equal(domain3D.facetFlags, flags_f3D)
        npt.assert_equal(domain2DRANS.vertexFlags, flags_v2DRANS)
        npt.assert_equal(domain2DRANS.segmentFlags, flags_s2DRANS)
        npt.assert_equal(domain3DRANS.vertexFlags, flags_v3DRANS)
        npt.assert_equal(domain3DRANS.facetFlags, flags_f3DRANS)


if __name__ == '__main__':

    unittest.main(verbosity=2)
