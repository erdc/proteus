SpatialTools
************

:doc:`BoundaryConditions`

Usage
=====

..
   General
   -------

Using spatial tools

#. Create a domain from `proteus.Domain`.
#. Create geometries that will be part of this domain from
   `proteus.SpatialTools` or `proteus.mprans.SpatialTools`.
#. Assemble domain.

Import
------

The base `SpatialTools` module can be imported as follows:

.. code-block:: python

   from proteus import SpatialTools as st

The `mprans` version of `SpatialTools`, adding functionality specific to
multi-phase flows, can be imported as follows:

.. code-block:: python

   from proteus.mprans import SpatialTools as st

Choosing the domain
-------------------

The domain class instance is what will hold all geometrical information once it
is assembled, and must be passed as an argument to all shapes created from
SpatialTools. Usually, the following class types are used from
`proteus.Domain`:

* 2D: ``PlanarStraightLineGraphDomain``
* 3D: ``PiecewiseLinearComplexDomain``

Those classes should be instantiated with no argument, e.g.:

.. code-block:: python

   from proteus import Domain
   domain = Domain.PlanarStraightLineGraphDomain().

Assembling the domain
---------------------

A very important final step not to forget is to assemble the domain once
all the desired geometries have been defined, so that the ``Domain`` instance
from `proteus.Domain` holds all the relevant geometrical information necessary
for running the simulation.

.. code-block:: python

   st.assembleDomain(domain)

BoundaryConditions with SpatialTools
------------------------------------

Usage
^^^^^

The mprans version of SpatialTools also associate geometry boundaries with
`BoundaryCondition` instances from mprans.BoundaryConditions. These
`BoundaryCondition` instances are easily accessible through a dictionary with
predefined names, e.g. for a `Rectangle` shape instance named `my_rectangle`,
``my_rectangle.BC['x-']`` is for accessing the boundary conditions of the left
segment, ``'x+'`` for the right, ``'y+'`` for the top, and ``'y-'`` for the
bottom.

Refer to :doc:`BoundaryConditions <./boundary_conditions>` documentation for
more information about the BoundaryConditions module.

.. note::

   The boundary conditions associated with the geometry do not have to be
   modified/set before assembling the domain, but removing/adding boundary
   conditions to a geometry on top of the predefined ones must be done before.

Linking it to _p.py files
^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   This does not apply to the `TwoPhaseFlow` module, which takes care of
   setting this automatically.

Linking the boundary conditions to the physical options _p.py files is done the following way (here for `RANS2P` boundary condition dictionnaries):

.. code-block:: python

   dirichletConditions = {0: lambda x, flag: domain.bc[flag].p_dirichlet.uOfXT,
                          1: lambda x, flag: domain.bc[flag].u_dirichlet.uOfXT,
                          2: lambda x, flag: domain.bc[flag].v_dirichlet.uOfXT,
                          3: lambda x, flag: domain.bc[flag].w_dirichlet.uOfXT}

   advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].p_advective.uOfXT,
                                      1: lambda x, flag: domain.bc[flag].u_advective.uOfXT,
                                      2: lambda x, flag: domain.bc[flag].v_advective.uOfXT,
                                      2: lambda x, flag: domain.bc[flag].w_advective.uOfXT}

   diffusiveFluxBoundaryConditions = {0:{},
                                      1:{1: lambda x, flag: domain.bc[flag].u_diffusive.uOfXT},
                                      2:{2: lambda x, flag: domain.bc[flag].v_diffusive.uOfXT},
                                      3:{3: lambda x, flag: domain.bc[flag].w_diffusive.uOfXT}}

This is always the same in the _p files, as long as it is pointing to the right
boundary conditions (e.g. p_dirichlet for pressure dirichlet). The boundary
conditions themselves can and should be manipulated externally (not from the
_p.py file), such as in the file where the geometries are first defined.


Complete Examples
=================

2D
---

.. code-block:: python

   from proteus import Domain
   from proteus.mprans import SpatialTools as st

   domain = Domain.PlanarStraightLineGraphDomain()

   my_tank = st.Tank2D(domain=domain,
                       dim=[10.,5.])

   my_rectangle = st.Rectangle(domain=domain,
                               dim=[1.,1.]
                               coords=[5.,2.5],
                               barycenter=[5.,2.5])
   my_rectangle.rotate(rot=3.14/4)
   my_rectangle.translate(trans=[0.1,0.1])

   st.assembleDomain(domain)

   my_tank.BC['x-'].setNoSlip()
   my_tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 0.1*x
   my_tank.BC['x-'].p_dirichlet.uOfXT = lambda x, t: -0.1*x
   my_tank.BC['x+'].setFreeSlip()
   my_tank.BC['y-'].setFreeSlip()
   my_tank.BC['y+'].setAtmosphere()

   my_rectangle.BC['x-'].setNoSlip()
   my_rectangle.BC['x+'].setNoSlip()
   my_rectangle.BC['y-'].setNoSlip()
   my_rectangle.BC['y+'].setNoSlip()
   

3D
---

.. code-block:: python

   from proteus import Domain
   from proteus.mprans import SpatialTools as st

   domain = Domain.PiecewiseLinearComplexDomain()

   my_tank = st.Tank3D(domain=domain,
                       dim=[10.,10.,5.])

   my_cylinder = st.Cylinder(domain=domain,
                             radius=1.,
                             height=3.,
                             nPoints=20,
                             coords=[5.,5.,2.5],
                             barycenter=[5.,5.,2.5])
   my_cylinder.rotate(rot=3.14/4,
                      axis=[1.,0.,0.],
                      pivot=my_cylinder.barycenter)
   my_cylinder.translate(trans=[0.1,0.1,0.1])

   st.assembleDomain(domain)

   my_tank.BC['x-'].setNoSlip()
   my_tank.BC['x-'].u_dirichlet.uOfXT = lambda x, t: 0.1*x
   my_tank.BC['x-'].p_dirichlet.uOfXT = lambda x, t: -0.1*x
   my_tank.BC['x+'].setFreeSlip()
   my_tank.BC['y-'].setFreeSlip()
   my_tank.BC['y+'].setFreeSlip()
   my_tank.BC['z-'].setFreeSlip()
   my_tank.BC['z+'].setFreeSlip()

   my_cylinder.BC['x-'].setNoSlip()
   my_cylinder.BC['x+'].setNoSlip()
   my_cylinder.BC['y-'].setNoSlip()
   my_cylinder.BC['y+'].setNoSlip()

Classes
=======


Base classes
------------

The following classes are accessible with an import from `proteus.SpatialTools`
and/or `proteus.mprans.SpatialTools`. Importing them from the mprans module
adds functionality such as the possibility to set multi-phase flow boundary
conditions and relaxation zones.

This is the same procedure as creating a `Domain` from scratch, with the added
benefit of being able to add more shapes to the domain as separate instances
and easy access to boundary conditions.

CustomShape
^^^^^^^^^^^

The most flexible type of shape, where everything is defined by the user. Any
geometry can be created with this. The minimum arguments necessary for setting
a custom geometry in 2D are: ``domain``, ``boundaryTags``, ``vertices``,
``vertexFlags``, ``segments``, and ``segmentFlags``. In 3D, the necessary
arguments are: ``domain``, ``boundaryTags``, ``vertices``, ``vertexFlags``,
``facets``, and ``facetFlags``. For additional arguments, please refer to the
source code in `proteus.SpatialTools`.

.. code-block:: python

   boundaryTags = {'my_tag1': 1,
                   'my_tag2': 2,
                   'my_tag3': 3}
   vertices = [[0.,0.],
               [1.,0.],
               [1.,1.],
               [0.,1.]]
   vertexFlags = [boundaryTags['my_tag1'],
                  boundaryTags['my_tag1'],
                  boundaryTags['my_tag2'],
                  boundaryTags['my_tag2']]
   segments = [[0, 1],
               [1, 2],
               [2, 3],
               [3, 0]]
   # flags can also be set from numbers included in the boundaryTags dictionary
   segmentFlags = [1, 2, 3, 2]
   my_customshape = st.CustomShape(domain=domain,
                                   vertices=vertices,
                                   vertexFlags=vertexFlags,
                                   segments=segments,
                                   segmentFlags=segmentFlags,
                                   boundaryTags=boundaryTags)
   my_customshape.BC['my_tag1'].setNoSlip()

Rectangle
^^^^^^^^^

A simple rectangular shape.

.. code-block:: python

    my_rectangle = st.Rectangle(domain=domain,
                                dim=[10.,2.],
                                coords=[5.,1.],
                                barycenter=[5.,1.])

Circle
^^^^^^

A simple circular shape.

.. code-block:: python

    my_circle = st.Circle(domain=domain,
                          radius=5.,
                          coords=[5.,5.],
                          barycenter=[5.,5.],
                          nPoints=20)

Cuboid
^^^^^^

A simple cuboidal shape.

.. code-block:: python

    my_cuboid = st.Cuboid(domain=domain,
                          dim=[10.,10.,2.],
                          coords=[5.,5.],
                          barycenter=[5.,5.])

Cylinder
^^^^^^^^

A simple cylindrical shape.

.. code-block:: python

    my_cylinder = st.Circle(domain=domain,
                            radius=5.,
                            height=10.
                            nPoints=20,
                            coords=[5.,5.,7.5],
                            barycenter=[5.,5.,7.5])

Sphere
^^^^^^

A simple spherical shape.

.. code-block:: python

    my_sphere = st.Sphere(domain=domain,
                          radius=5.,
                          coords=[2.,2.],
                          barycenter=[2.,2.],
                          nSectors=10)

ShapeSTL
^^^^^^^^

For importing STL geometries. It needs a `.stl` ASCII file, and does not
currently work with binary files. The STL geometry is converted in a Proteus
readable format, automatically creating vertices and facets, and a single
boundary tag/flag for the whole STL geometry.

.. code-block:: python

    my_stl = st.ShapeSTL(domain=domain,
                         filename='path/to/my/file.stl')

mprans specific
---------------

The following classes are for use with multi-phase flow and can only be
imported from `proteus.mprans.SpatialTools`.


Tank2D
^^^^^^

The `Tank2D` class can be used to create a rectangular tank. This class allows
for "sponge layers", or "relaxation zones" that are usually used for wave
absorption or wave generation to get rid of reflected waves in the domain. The
lower left corner of the tank is at the origin `[0.,0.]` when created (but it
can still be translated later on), and sponge layers extend outwards of the
numerical tank. A `Tank2D` of dimensions [10.,2.] and sponge layers of
length 3. on both sides will have a total domain size of [16,2], spanning
from x=-3 to x=13.

.. code-block:: python

   my_tank = st.Tank2D(domain=domain,
                       dim=[10.,2.])
   # make sponge layers
   my_tank.setSponge(x_n=3., x_p=3.)
   # set absorption zone (x_p -> x+)
   my_tank.setAbsorptionZones(dragAlpha=1.e6,
                              x_p=True)
   # set generation zone
   from proteus import WaveTools as wt (x_n -> x-)
   my_wave = wt.MonochromaticWave()
   he = 0.01
   my_tank.setGenerationZones(dragAlpha=1.e6,
                              smoothing=3*he,
                              wave=wave,
                              x_n=True)
   # set boundary conditions
   my_tank.BC['y+'].setAtmosphere()
   my_tank.BC['y-'].setFreeSlip()
   my_tank.BC['x+'].setFreeSlip()
   my_tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=my_wave
                                                     smoothing=3*he)
   my_tank.BC['sponge'].setNonMaterial()

.. important::

   `Tank2D` instances should not be rotated as this can lead to problems with
   relaxation zones and boundary conditions.


Tank3D
^^^^^^

Very similar to the `Tank3D`, it is a cuboid for 3D domains with the
possibility of adding sponge layers.

.. code-block:: python

   my_tank = st.Tank2D(domain=domain,
                       dim=[10.,10., 2.])
   # make sponge layers
   my_tank.setSponge(x_n=3., x_p=3., y_p=3., y_n=3.)
   # set absorption zones
   my_tank.setAbsorptionZones(dragAlpha=1.e6,
                              x_p=True,
                              y_p=True,
                              y_n=True)
   # set generation zone
   from proteus import WaveTools as wt (x_n -> x-)
   my_wave = wt.MonochromaticWave()
   he = 0.01
   my_tank.setGenerationZones(dragAlpha=1.e6,
                              smoothing=3*he,
                              wave=wave,
                              x_n=True)
   # set boundary conditions
   my_tank.BC['z+'].setAtmosphere()
   my_tank.BC['z-'].setFreeSlip()
   my_tank.BC['y-'].setFreeSlip()
   my_tank.BC['x+'].setFreeSlip()
   my_tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave=my_wave
                                                     smoothing=3*he)
   my_tank.BC['sponge'].setNonMaterial()

.. important::

   `Tank3D` instances should not be rotated as this can lead to problems with
   relaxation zones and boundary conditions.


TankWithObstacle2D
^^^^^^^^^^^^^^^^^^
