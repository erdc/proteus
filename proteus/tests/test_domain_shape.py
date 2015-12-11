"""
Testing module for Domain.py, Shape.py, BC.py
Work in progress
"""

def checkFlags(domain):
    """
    Function to check that local flags match domain flags
    :arg domain: domain (with list of shapes)
    """
    flag_nb = 0
    i_v = 0
    i_s = 0
    i_f = 0
    i_bc = 1  # the first bc (index 0) in domain if ignored
    assert len(domain.shape_list) != 0, \
        'Domain does not have Shapes associated to it!'
    for shape in domain.shape_list:
        for flag in shape.vertexFlags:
            assert flag+flag_nb == domain.vertexFlags[i_v], \
                'vertex flags not matching!'
            i_v += 1
        if shape.segments is not None:
            for flag in shape.segmentFlags:
                assert flag+flag_nb == domain.segmentFlags[i_s], \
                    'segment flags not matching!'
                i_s += 1
        if shape.facets is not None:
            for flag in shape.facetFlags:
                assert flag+flag_nb == domain.facetFlags[i_f], \
                    'facet flags not matching!'
                i_f += 1
        for bc in shape.BC_list:
            assert bc == domain.bc[i_bc], \
                    'boundary conditions not matching!'
            i_bc += 1
        flag_nb += len(shape.BC_list)
    nd = str(int(domain.nd))
    name = domain.name
    print(name + ' and its ' + nd + 'D shapes passed flag indexing test')



if __name__ == '__main__':

    print('|-----------------------------------------|')
    print('| Testing SpatialTools and Domain modules |')
    print('|-----------------------------------------|')

    from proteus import Domain
    from proteus import SpatialTools as st
    from proteus.default_n import *
    # ----- DOMAIN ----- #
    domain3D = Domain.PiecewiseLinearComplexDomain()
    domain3D.Mesh.elementSize(1.)
    domain2D = Domain.PlanarStraightLineGraphDomain()
    domain2D.Mesh.elementSize(1.)
    # ----- SHAPES ----- #
    tank3D = st.Tank3D(domain3D, [10., 10., 10.])
    caisson3D = st.Cuboid(domain3D, [1., 1., 1.], [2., 2., 2.])
    rigidbody3D = st.BodyCuboid(domain3D, [1., 1., 1.], [4., 4., 4.])
    tank2D = st.Tank2D(domain2D, [10., 10.])
    caisson2D = st.Rectangle(domain2D, [1., 1.], [2., 2.])
    rigidbody2D = st.BodyRectangle(domain2D, [1., 1.], [4., 4.])
    # ----- CUSTOM SHAPE ----- #
    bt3D = {'bottom': 1, 'front':2, 'right':3, 'back': 4, 'left':5,
            'top':6, 'obstacle':7}
    vertices3D = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
                  [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]]
    vertexFlags3D = [bt3D['left'], bt3D['right'], bt3D['right'], bt3D['left'],
                     bt3D['left'], bt3D['right'], bt3D['right'], bt3D['left']]
    facets3D = [[[0, 1, 2, 3]], [[0, 1, 5, 4]], [[1, 2, 6, 5]], [[2, 3, 7, 6]],
                [[3, 0, 4, 7]], [[4, 5, 6, 7]]]
    facetFlags3D = [bt3D['bottom'], bt3D['front'], bt3D['right'], bt3D['back'],
                    bt3D['left'], bt3D['top']]
    custom3D = st.CustomShape(domain=domain3D,
                              vertices=vertices3D,
                              vertexFlags=vertexFlags3D,
                              facets=facets3D,
                              facetFlags=facetFlags3D,
                              boundaryTags=bt3D)
    bt2D = {'bottom': 5, 'right': 4, 'left': 7, 'top': 6}
    vertices2D = [[0., 0.], [1., 0.], [1., 1.], [0., 1.]]
    vertexFlags2D = [bt2D['bottom'], bt2D['bottom'], bt2D['top'], bt2D['top']]
    segments2D = [[0, 1], [1, 2], [2, 3], [3, 0]]
    segmentFlags2D = [bt2D['bottom'], bt2D['right'], bt2D['top'], bt2D['left']]
    custom3D = st.CustomShape(domain=domain2D,
                              vertices=vertices2D,
                              vertexFlags=vertexFlags2D,
                              facets=segments2D,
                              facetFlags=segmentFlags2D,
                              boundaryTags=bt2D)

    # ----- BUILD DOMAINS ----- #
    st.buildDomain(domain3D)
    st.buildDomain(domain2D)
    # --------- TESTS --------- #
    checkFlags(domain3D)
    checkFlags(domain2D)
