"""
Module creating predifined shapes associated. Each shape needs a domain as argument and the wanted dimensions.
Boundary conditions objects are automatically created for each facet (3D) or segment (2D) defining the shape.
classes:
- Shape: super class, regroups functions common to all shapes
- Cuboid: creates a cuboid 
- Rectangle: creates a rectangle

TO DO:

Add more predifined shape types
"""

import BC as bc
import numpy as np

class Shape:
    """
    Class defining a shape
    :param domain: domain in which the shape is defined
    """
    def __init__(self, domain):
        self.domain = domain
        self.vertices = None
        self.vertexFlags = None
        self.segments = None
        self.segmentFlags = None
        self.facets = None
        self.facetFlags = None
        self.regions = None
        self.regionFlags = None
        self.holes = None
        self.volume = None
        self.mass = None
        self.density = None
        self.barycenter = None
        self.cm = None  # centre of mass
        self.bc = None
        self.snv = len(self.domain.vertices)  # total number of vertices in domain when shape.__init__
        self.snf = len(self.domain.facets)
        self.sns = len(self.domain.segments)
        self.snr = len(self.domain.regions)
        self.snbc = len(self.domain.bc)

    def addShape(self):
        # add new information to the Shape class
        self.domain.vertices += self.vertices  # append new vertices to static variable of Shape
        self.domain.vertexFlags += [vertexFlag+self.snv for vertexFlag in self.vertexFlags]
        self.domain.segments += [[vertex_nb+self.snv for vertex_nb in segment] for segment in self.segments]
        self.domain.segmentFlags += [segmentFlag+self.snv for segmentFlag in self.segmentFlags]
        if self.domain.nd == 3:
            self.domain.facets += [[[vertex_nb+self.snv for vertex_nb in facet]] for facet in self.facets]
            self.domain.facetFlags += [facetFlag+self.snf for facetFlag in self.facetFlags]
            self.domain.holes += self.holes
        self.domain.regions += self.regions
        self.domain.regionFlags += [regionFlag+self.snf for regionFlag in self.regionFlags]
        self.domain.bc += self.bc
        self.domain.update()

    def setMass(self, mass):
        self.mass = mass
        if self.volume is not None:
            self.density = self.mass/self.volume

    def setDensity(self, density):
        self.density = density
        if self.volume is not None:
            self.mass = self.density*self.volume

    def setPosition(self, coords):
        self.coords = coords
        vertices = np.asarray(self.vertices)
        regions = np.asarray(self.regions)
        vertices = vertices - self.coords + coords
        regions = regions - self.regions + coords
        self.vertices = vertices.tolist()
        self.regions = regions.tolist()
        self.updateDomain()

    def updateDomain(self):
        nv = len(self.vertices)
        nr = len(self.regions)
        self.domain.vertices[self.snv:self.snv+nv] = self.vertices
        self.domain.regions[self.snr:self.snr+nr] = self.regions
        self.domain.update()

    def setConstraints(self, free_x, free_r):
        self.free_x = free_x
        self.free_r = free_r


class Cuboid(Shape):
    def __init__(self, domain, dim, coords=(0.,0.,0.)):
        Shape.__init__(self, domain)
        self.dim = dim  # length, width height
        self.volume = dim[0]*dim[1]*dim[2]
        self.coords = coords
        self.coords = x, y, z = coords
        self.vertices = [[x-0.5*dim[0], y-0.5*dim[1], z-0.5*dim[2]],
                         [x-0.5*dim[0], y+0.5*dim[1], z-0.5*dim[2]],
                         [x+0.5*dim[0], y+0.5*dim[1], z-0.5*dim[2]],
                         [x+0.5*dim[0], y-0.5*dim[1], z-0.5*dim[2]],
                         [x-0.5*dim[0], y-0.5*dim[1], z+0.5*dim[2]],
                         [x-0.5*dim[0], y+0.5*dim[1], z+0.5*dim[2]],
                         [x+0.5*dim[0], y+0.5*dim[1], z+0.5*dim[2]],
                         [x+0.5*dim[0], y-0.5*dim[1], z+0.5*dim[2]]]
        self.segments = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6],
                         [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]]
        self.facets = [[0, 1, 2, 3],  # bottom
                       [0, 1, 5, 4],  # front
                       [1, 2, 6, 5],  # right
                       [2, 3, 7, 6],  # back
                       [3, 0, 4, 7],  # left
                       [4, 5, 6, 7]]  # top
        self.regions = [[x, y, z]]
        self.holes = [coords]
        # defining flags for boundary conditions
        self.facetFlags = [0, 1, 2, 3, 4, 5]  # bottom, front, right, back, left, top
        self.vertexFlags = [0, 0, 0, 0, 5, 5, 5, 5]  # only top and bottom for vertices
        self.segmentFlags = []
        self.regionFlags = [0]
        # Initialize (empty) boundary conditions
        self.bc_bottom = bc.BoundaryConditions()
        self.bc_front = bc.BoundaryConditions()
        self.bc_right = bc.BoundaryConditions()
        self.bc_back = bc.BoundaryConditions()
        self.bc_left = bc.BoundaryConditions()
        self.bc_top = bc.BoundaryConditions()
        self.bc = [self.bc_bottom, self.bc_front, self.bc_right,
                    self.bc_back, self.bc_left, self.bc_top]
        self.addShape()  # adding shape to domain


class Rectangle(Shape):
    def __init__(self, domain, dim, coords=(0.,0.)):
        Shape.__init__(self, domain)
        self.dim = dim
        self.coords = x, y = coords
        self.vertices = [[x-0.5*dim[0], y-0.5*dim[1]],
                         [x+0.5*dim[0], y-0.5*dim[1]],
                         [x+0.5*dim[0], y+0.5*dim[1]],
                         [x-0.5*dim[0], y+0.5*dim[1]]]
        self.segments= [[0,1],[1,2],[2,3],[3,0]]
        self.regions = [[x, y]]
        self.segmentFlags = [0, 1, 2, 3] # bottom, right, top, left
        self.vertexFlags = [0, 0, 2, 2]  # bottom, bottom, top, top
        self.regionFlags = [0]
        # Initialize (empty) boundary conditions
        self.bc_bottom = bc.BoundaryConditions()
        self.bc_right = bc.BoundaryConditions()
        self.bc_left = bc.BoundaryConditions()
        self.bc_top = bc.BoundaryConditions()
        self.bc = [self.bc_bottom, self.bc_right, self.bc_left, self.bc_top]
        self.addShapes()  # adding shape to domain
