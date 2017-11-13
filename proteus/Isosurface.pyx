#!python
#cython: embedsignature=True
"""
AuxiliaryVariables subclasses for extracting isosurfaces and contours
"""
from collections import defaultdict, OrderedDict
from itertools import product
import os
#from proteus.EGeometry import etriple, ecross, enorm, edot

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
cimport numpy as np
from numpy.linalg import norm

from . import Comm
from .AuxiliaryVariables import AV_base
from proteus import Profiling
from .Profiling import logEvent as log

from libc.math cimport sqrt

DEF X=0
DEF Y=1
DEF Z=2

cdef np.ndarray[np.float64_t,ndim=1] EVec(double x=0.0, double y=0.0, double z=0.0):
    cdef  np.ndarray [np.float64_t, ndim=1] v = np.zeros((3,),'d')
    v[X] = x
    v[Y] = y
    v[Z] = z
    return v

cdef double enorm(np.ndarray[np.float64_t, ndim=1] v):
    return sqrt(v[X]**2 + v[Y]**2 + v[Z]**2)

cdef double edot(np.ndarray[np.float64_t, ndim=1] v0, np.ndarray[np.float64_t, ndim=1] v1):
    return v0[X]*v1[X] + v0[Y]*v1[Y] + v0[Z]*v1[Z]

cdef np.ndarray[np.float64_t, ndim=1] ecross(np.ndarray[np.float64_t, ndim=1] v0,  np.ndarray[np.float64_t, ndim=1] v1):
    return EVec(v0[Y]*v1[Z] - v0[Z]*v1[Y],
                v0[Z]*v1[X] - v0[X]*v1[Z],
                v0[X]*v1[Y] - v0[Y]*v1[X])

cdef double etriple(np.ndarray[np.float64_t, ndim=1] v0, np.ndarray[np.float64_t, ndim=1] v1, np.ndarray[np.float64_t, ndim=1] v2):
    return edot(v0,ecross(v1,v2))


class Isosurface(AV_base):

    """Extract isosurfaces"""

    def __init__(self, isosurfaces, domain, activeTime=None,
                 sampleRate=0, format='pov', writeBoundary=True):
        """Create a set of isosurfaces that will be extracted and serialized

        :param isosurfaces: An iterable of "isosurfaces".  Each isosurface is
        specified by a 2-tuple, with the first element in the tuple is a field
        from which to extract isosurfaces, and the second element is an n-tuple
        of isosurface values to extract.

        :param domain: a Domain object

        :param activeTime: If not None, a 2-tuple of start time and end time
        for which the point gauge is active.

        :param sampleRate: The intervals at which samples should be measured.
        Note that this is a rough lower bound, and that the gauge values could
        be computed less frequently depending on the time integrator.  The
        default value of zero computes the gauge values at every time step.

        :param format: the file format for the isosurfaces ('pov', None).

        :param writeBoundary: whether to write the boundary mesh to or the
        domain bounding box

        Example:

        phi0 = (isosurfaces=('phi', (0.0,)),
                activeTime=(0, 2.5),
                sampleRate=0.2,
                fileName='water.pov')

        This creates an Isosurfaces object that will extract phi(x,y,z,t) = 0
        at simulation time between = 0 and 2.5 with samples taken no more
        frequently than every 0.2 seconds.  Results will be saved to:
        water.pov.

        """
        self.isosurfaces = isosurfaces
        self.domain = domain
        self.nodes = {}
        self.elements = {}
        self.normals = {}
        self.normal_indices = {}
        for f, values in self.isosurfaces:
            for v in values:
                assert v == 0.0, \
                    "only implemented for 0 isosurface in 3D for now"
        self.activeTime = activeTime
        self.sampleRate = sampleRate
        self.format = format
        self.comm = Comm.get()
        self.writeBoundary = writeBoundary
        self.fileprefix = 'isosurface'

    def attachModel(self, model, ar):
        """ Attach this isosurface to the given simulation model.
        """
        self.model = model
        fine_grid = model.levelModelList[-1]
        self.fieldNames = fine_grid.coefficients.variableNames
        self.elementNodesArray = fine_grid.mesh.elementNodesArray
        self.exteriorElementBoundariesArray \
            = fine_grid.mesh.exteriorElementBoundariesArray
        self.elementBoundaryNodesArray \
            = fine_grid.mesh.elementBoundaryNodesArray
        self.elementBoundaryMaterialTypes \
            = fine_grid.mesh.elementBoundaryMaterialTypes
        self.boundaryNodes = set(
            self.elementBoundaryNodesArray[
                self.exteriorElementBoundariesArray
            ].flatten())
        self.g2b = dict((gn, bn) for bn, gn in enumerate(self.boundaryNodes))
        self.nodeArray = fine_grid.mesh.nodeArray
        self.num_owned_elements = fine_grid.mesh.nElements_global
        self.u = fine_grid.u
        self.timeIntegration = fine_grid.timeIntegration
        self.nFrames = 0
        self.next_output = 0
        self.writeSceneHeader()
        return self

    def attachHDF5(self, h5, step, cam=None):
        """
        Attach this isosurface to and HDF5 archive
        """
        from collections import namedtuple
        self.fieldNames = [isosurface[0] for isosurface in self.isosurfaces]
        print("ATTACHING TO HDF5 !!", self.fieldNames)
        self.elementNodesArray = h5.get_node("/elementsSpatial_Domain" +
                                            repr(step))[:]
        self.nodeArray = h5.get_node("/nodesSpatial_Domain" + repr(step))[:]
        self.num_owned_elements = len(self.elementNodesArray)
        self.u = {}
        FemField = namedtuple('FemField', ['dof'])
        for field_i, field in enumerate(self.fieldNames):
            self.u[field_i] = FemField(dof=h5.get_node("/" +
                                                      self.isosurfaces[0][0] +
                                                      repr(step))[:])
        self.nFrames = step
        self.next_output = 0
        if step == 0:
            self.writeSceneHeader(cam)
        return self

    def triangulateIsosurface(self, field, value):
        """
        Build a triangular mesh of the isosurface
        """
        cdef int eN, i, J, I, nMinus, nPlus, nZeros, nN_start
        cdef double eps, s
        cdef np.ndarray[np.float64_t, ndim=1] x, normal
        self.nodes[(field, value)] = nodes = []
        self.elements[(field, value)] = elements = []
        self.normals[(field, value)] = normals = []
        self.normal_indices[(field, value)] = normal_indices = []
        if value != 0.0:
            raise NotImplementedError("Only zero isocontour extraction")
        phi = self.u[self.fieldNames.index(field)].dof
        for eN in range(self.num_owned_elements):
            plus = []
            minus = []
            zeros = []
            eps = 1.0e-8
            for i in range(4):
                I = self.elementNodesArray[eN, i]
                if phi[I] > eps:
                    plus.append(I)
                elif phi[I] < -eps:
                    minus.append(I)
                else:
                    zeros.append(I)
            nZeros = len(zeros)
            nMinus = len(minus)
            nPlus = len(plus)
            assert(nZeros + nMinus + nPlus == 4)
            nN_start = len(nodes)
            if nMinus == 2 and nPlus == 2:  # 4 cut edges
                for J in minus:
                    for I in plus:
                        s = -phi[I] / (phi[J] - phi[I])
                        x = s * (self.nodeArray[J] - self.nodeArray[I]) \
                            + self.nodeArray[I]
                        nodes.append(x)
                elements.append([nN_start + j for j in range(3)])
                elements.append([nN_start + j + 1 for j in range(3)])
                if etriple(self.nodeArray[plus[0]] - nodes[-4],
                           nodes[-3] - nodes[-4],
                           nodes[-2] - nodes[-4]) < 0.0:
                    elements[-2] = [elements[-2][0],
                                    elements[-2][2],
                                    elements[-2][1]]
                if etriple(self.nodeArray[plus[1]] - nodes[-3],
                           nodes[-2] - nodes[-3],
                           nodes[-1] - nodes[-3]) < 0.0:
                    elements[-1] = [elements[-1][0],
                                    elements[-1][2],
                                    elements[-1][1]]
                normal = ecross(
                    nodes[elements[-2][1]] - nodes[elements[-2][0]],
                    nodes[elements[-2][2]] - nodes[elements[-2][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-2])
                normal = ecross(
                    nodes[elements[-1][1]] - nodes[elements[-1][0]],
                    nodes[elements[-1][2]] - nodes[elements[-1][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-1])
            elif nPlus == 3 and nMinus == 1:  # 3 cut edges
                I = minus[0]
                for J in plus:
                    s = -phi[I] / (phi[J] - phi[I])
                    x = s * (self.nodeArray[J] - self.nodeArray[I]) + \
                        self.nodeArray[I]
                    nodes.append(x)
                elements.append([nN_start + j for j in range(3)])
                if etriple(self.nodeArray[minus[0]] - nodes[-3],
                           nodes[-2] - nodes[-3],
                           nodes[-1] - nodes[-3]) > 0.0:
                    elements[-1] = [elements[-1][0],
                                    elements[-1][2],
                                    elements[-1][1]]
                normal = ecross(
                    nodes[elements[-1][1]] - nodes[elements[-1][0]],
                    nodes[elements[-1][2]] - nodes[elements[-1][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-1])
            elif nMinus == 3 and nPlus == 1:  # 3 cut edges
                I = plus[0]
                for J in minus:
                    s = -phi[I] / (phi[J] - phi[I])
                    x = s * (self.nodeArray[J] - self.nodeArray[I]) + \
                        self.nodeArray[I]
                    nodes.append(x)
                elements.append([nN_start + j for j in range(3)])
                if etriple(self.nodeArray[plus[0]] - nodes[-3],
                           nodes[-2] - nodes[-3],
                           nodes[-1] - nodes[-3]) < 0.0:
                    elements[-1] = [elements[-1][0],
                                    elements[-1][2],
                                    elements[-1][1]]
                normal = ecross(
                    nodes[elements[-1][1]] - nodes[elements[-1][0]],
                    nodes[elements[-1][2]] - nodes[elements[-1][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-1])
            elif nZeros == 1 and ((nPlus == 2 and nMinus == 1) or
                                  (nPlus == 1 and nMinus == 2)):
                # 2 cut edges, 1 vertex lies in plane
                nodes.append(self.nodeArray[zeros[0]])
                for J in minus:
                    for I in plus:
                        s = -phi[I] / (phi[J] - phi[I])
                        x = s * (self.nodeArray[J] - self.nodeArray[I]) \
                            + self.nodeArray[I]
                        nodes.append(x)
                elements.append([nN_start + j for j in range(3)])
                if etriple(self.nodeArray[plus[0]] - nodes[-3],
                           nodes[-2] - nodes[-3],
                           nodes[-1] - nodes[-3]) < 0.0:
                    elements[-1] = [elements[-1][0],
                                    elements[-1][2],
                                    elements[-1][1]]
                normal = ecross(
                    nodes[elements[-1][1]] - nodes[elements[-1][0]],
                    nodes[elements[-1][2]] - nodes[elements[-1][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-1])
            elif nZeros == 2 and nPlus == 1 and nMinus == 1:
                # 1 cut edge, 2 vertices lie in plane
                for I in zeros:
                    nodes.append(self.nodeArray[I])
                I = plus[0]
                J = minus[0]
                s = -phi[I] / (phi[J] - phi[I])
                x = s * (self.nodeArray[J] - self.nodeArray[I]) \
                    + self.nodeArray[I]
                nodes.append(x)
                elements.append([nN_start + j for j in range(3)])
                if etriple(self.nodeArray[plus[0]] - nodes[-3],
                           nodes[-2] - nodes[-3],
                           nodes[-1] - nodes[-3]) < 0.0:
                    elements[-1] = [elements[-1][0],
                                    elements[-1][2],
                                    elements[-1][1]]
                normal = ecross(
                    nodes[elements[-1][1]] - nodes[elements[-1][0]],
                    nodes[elements[-1][2]] - nodes[elements[-1][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-1])
            elif nZeros == 3:  # 3 vertices lie in plane
                for I in zeros:
                    nodes.append(self.nodeArray[I])
                elements.append([nN_start + j for j in range(3)])
                if nPlus == 1:
                    if etriple(self.nodeArray[plus[0]] - nodes[-3],
                               nodes[-2] - nodes[-3],
                               nodes[-1] - nodes[-3]) < 0.0:
                        elements[-1] = [elements[-1][0],
                                        elements[-1][2],
                                        elements[-1][1]]
                else:
                    if etriple(self.nodeArray[minus[0]] - nodes[-3],
                               nodes[-2] - nodes[-3],
                               nodes[-1] - nodes[-3]) > 0.0:
                        elements[-1] = [elements[-1][0],
                                        elements[-1][2],
                                        elements[-1][1]]
                normal = ecross(
                    nodes[elements[-1][1]] - nodes[elements[-1][0]],
                    nodes[elements[-1][2]] - nodes[elements[-1][0]])
                normal /= enorm(normal)
                normals.append(normal)
                normals.append(normal)
                normals.append(normal)
                normal_indices.append(elements[-1])

    def writeIsosurfaceMesh(self, field, value, frame):
        if self.format == 'pov':
            log("Writing pov frame " + repr(frame))
            self.writeIsosurfaceMesh_povray(field, value, frame)
        elif self.format == 'h5':
            self.writeIsosurfaceMesh_h5(field, value, frame)
        elif self.format is None:
            pass
        else:
            log("Isosurface file format not recognized")

    def writeIsosurfaceMesh_h5(self, field, value, frame):
        import h5py
        nodes = self.nodes[(field, value)]
        elements = self.elements[(field, value)]
        normals = self.normals[(field, value)]
        normal_indices = self.normal_indices[(field, value)]
        filename = os.path.join(Profiling.logDir, self.fileprefix+str(self.comm.rank())+'.h5')
        if self.nFrames == 0:
            f = h5py.File(filename, "w")
        else:
            f = h5py.File(filename, "a")
        dset = f.create_dataset('nodes'+str(self.nFrames), data=nodes)
        dset = f.create_dataset('elems'+str(self.nFrames), data=elements)
        dset = f.create_dataset('normals'+str(self.nFrames), data=normals)
        dset = f.create_dataset('normal_indices'+str(self.nFrames), data=normal_indices)
        f.close()

    def writeIsosurfaceMesh_povray(self, field, value, frame):
        """
        Write the triangular mesh to a povray file
        """
        from string import Template
        nodes = self.nodes[(field, value)]
        elements = self.elements[(field, value)]
        normals = self.normals[(field, value)]
        normal_indices = self.normal_indices[(field, value)]
        pov_filename = "{field:s}_{value:f}_{frame:04d}.pov".format(
            field=field,
            value=value,
            frame=self.nFrames)
        self.comm.beginSequential()
        if self.comm.isMaster():
            pov = open(pov_filename, "w")
            dx = np.array(self.domain.L) * 0.02
            nll = np.array(self.domain.x) - dx
            fur = np.array(self.domain.L) + dx
            fur[2] = self.domain.L[2] - dx[2]  # clip  the  top off
            pov.write("""#version 3.7;
#include  "proteus.inc"
""")
            if not self.writeBoundary:
                pov.write(Template("""
object
{
difference
{
box {
    <$nll_x,$nll_y,$nll_z>,  // Near lower left corner
    <$fur_x,$fur_y,$fur_z>   // Far upper right corner
    }
box {
    <$domain_nll_x,$domain_nll_y,$domain_nll_z>,  // Near lower left corner
    <$domain_fur_x,$domain_fur_y,$domain_fur_z>   // Far upper right corner
    }
}//difference of perturbed bounding box and boundary
                    matrix < 1.000000, 0.000000, 0.000000,
                             0.000000, 1.000000, 0.000000,
                             0.000000, 0.000000, 1.000000,
                             0.000000, 0.000000, 0.000000 >
                    tank_material()
}//object
""").substitute(nll_x=nll[0],
                    nll_y=nll[1],
                    nll_z=nll[2],
                    fur_x=fur[0],
                    fur_y=fur[1],
                    fur_z=fur[2],
                    domain_nll_x=self.domain.x[0],
                    domain_nll_y=self.domain.x[1],
                    domain_nll_z=self.domain.x[2],
                    domain_fur_x=self.domain.L[0],
                    domain_fur_y=self.domain.L[1],
                    domain_fur_z=self.domain.L[2]))
                pov.flush()
            if self.writeBoundary:
                pov.write(Template("""
object
{
//difference
//{
// box {
//    <$nll_x,$nll_y,$nll_z>,  // Near lower left corner
//    <$fur_x,$fur_y,$fur_z>   // Far upper right corner
//    }
union
{
""").substitute(nll_x=nll[0],
                    nll_y=nll[1],
                    nll_z=nll[2],
                    fur_x=fur[0],
                    fur_y=fur[1],
                    fur_z=fur[2]))
            pov.flush()
        else:
            pov = open(pov_filename, "a")
        if self.writeBoundary:
            povScene = """mesh2 {
vertex_vectors {"""
            pov.write(povScene)
            pov.write("{0:d},\n".format(len(self.boundaryNodes)))
            for n in self.boundaryNodes:
                pov.write("<{0:f}, {1:f}, {2:f}>,\n".format(
                    *self.nodeArray[n]))
            pov.write("}\n")
            pov.write("""face_indices {
                             """)
            pov.write("{0:d},\n".format(
                np.count_nonzero(self.elementBoundaryMaterialTypes)))
            for ebN in self.exteriorElementBoundariesArray:
                if self.elementBoundaryMaterialTypes[ebN] > 0:
                    bnt = tuple(
                        [self.g2b[n] for n
                         in self.elementBoundaryNodesArray[ebN]])
                    pov.write("<{0:d},{1:d},{2:d}>,\n".format(*bnt))
            pov.write("}\n")
            pov.write("//inside_vector on\n}//mesh\n")
            pov.flush()
            pov.close()
            self.comm.endSequential()
            self.comm.barrier()
            self.comm.beginSequential()
            pov = open(pov_filename, "a")
            if self.comm.isMaster():
                pov.write("""}//union of meshes
//}//difference of perturbed bounding box and boundary
                    matrix < 1.000000, 0.000000, 0.000000,
                             0.000000, 1.000000, 0.000000,
                             0.000000, 0.000000, 1.000000,
                             0.000000, 0.000000, 0.000000 >
                    tank_material()
}//object
""")
        povScene = """mesh2 {
vertex_vectors {"""
        pov.write(povScene)
        pov.write("{0:d},\n".format(len(nodes)))
        for n in nodes:
            pov.write("<{0:f}, {1:f}, {2:f}>,\n".format(*n))
        pov.write("""        }
                normal_vectors {
                        """)
        pov.write("{0:d},\n".format((len(normals))))
        for n in normals:
            pov.write("<{0:f}, {1:f}, {2:f}>,\n".format(*n))
        pov.write("""        }
                 face_indices {
                         """)
        pov.write("{0:d},\n".format(len(elements)))
        for e in elements:
            pov.write("<{0:d}, {1:d}, {2:d}>,\n".format(*e))
        pov.write("""        }
                 normal_indices {
                         """)
        pov.write("{0:d},\n".format(len(normal_indices)))
        for ni in normal_indices:
            pov.write("<{0:d}, {1:d}, {2:d}>,\n".format(*ni))
        pov.write("""        }
                    matrix < 1.000000, 0.000000, 0.000000,
                             0.000000, 1.000000, 0.000000,
                             0.000000, 0.000000, 1.000000,
                             0.000000, 0.000000, 0.000000 >
                    isosurface_material()
                }
""")
        pov.flush()
        pov.close()
        self.comm.endSequential()

    def calculate(self, checkTime=True):
        """Extracts isosourfaces at current time and update open output files
        """
        from string import Template
        if checkTime:
            time = self.timeIntegration.tLast
            log("Calculate called at time " + repr(time))
            # check that gauge is in its active time region
            if self.activeTime is not None and (self.activeTime[0] > time or
                                                self.activeTime[1] < time):
                return

            # check that gauge is ready to be sampled again
            if (time < self.next_output):
                return
        assert self.elementNodesArray.shape[1] == 4, \
            "Elements have {0:d} vertices but algorithm is for tets".format(
                self.elementNodesArray.shape[1])
        for field, values in self.isosurfaces:
            for v in values:
                self.triangulateIsosurface(field, v)
                self.writeIsosurfaceMesh(field, v, self.nFrames)
        self.nFrames += 1
        if checkTime and time != 0:
            self.next_output += self.sampleRate

    def writeSceneHeader(self, cam = None):
        """
        Write a scene description (can be modified before running povray)
        """
        from string import Template
        if self.comm.isMaster():
            look_at = [0.5 * (x + L)
                       for x, L in zip(self.domain.x, self.domain.L)]
            if cam is None:
                cam = [0.5*(self.domain.x[0] + self.domain.L[0]),
                       self.domain.x[1] - 2*self.domain.L[1],
                       self.domain.x[2] + 0.85*(self.domain.L[2] + self.domain.x[2])]
            light = [0.5 * (x + L)
                     for x, L in zip(self.domain.x, self.domain.L)]
            light[2] = self.domain.x[2] + 5 * self.domain.L[2]
            light_dx, light_dy = (self.domain.x[0] - light[0],
                                  self.domain.x[1] - light[1])
            floor_z = self.domain.x[2] - 0.01 * \
                self.domain.L[2]  # offset slightly
            wall_y = self.domain.x[1] - 2.0 * \
                self.domain.L[1]  # offset slightly
            sky_z = self.domain.x[2] + 10 * self.domain.L[2]
            pov = open("proteus.inc", "w")
            povSceneTemplate = Template("""#include "colors.inc"
#include "textures.inc"
#include "glass.inc"
#include "metals.inc"
#include "golds.inc"
#include "stones.inc"
#include "woods.inc"
#include "shapes.inc"
#include "shapes2.inc"
#include "functions.inc"
#include "math.inc"
#include "transforms.inc"

global_settings {
        ambient_light color rgb <1.0, 1.0, 1.0>
        assumed_gamma 2
}

background { color rgb <0.319997, 0.340002, 0.429999>}

camera {
        perspective
        location <$cam_x,$cam_y,$cam_z>
        sky <0.0, 0.0, 5.0>
        up <0, 0, 1>
        right <1.33, 0, 0>
        angle 45.000000
        look_at <$look_at_x,$look_at_y,$look_at_z>
}

light_source {<$light_x,$light_y,$light_z> color White}

light_source {
        <$light_x+$light_dx,$light_y+$light_dy,$light_z>
        color <0.99980, 0.999800, 0.999800>*2.250000
        spotlight
        point_at <0.5,0.5,0.0>
}

light_source {
        <$light_x+$light_dx,$light_y-$light_dy,$light_z>
        color <0.99980, 0.999800, 0.999800>*2.250000
        spotlight
        point_at <0.5,0.5,0.0>
}

light_source {
        <$light_x-$light_dx,$light_y-$light_dy,$light_z>
        color <0.99980, 0.999800, 0.999800>*2.250000
        spotlight
        point_at <0.5,0.5,0.0>
}

light_source {
        <$light_x-$light_dx,+$light_y,$light_z>
        color <0.99980, 0.999800, 0.999800>*2.250000
        spotlight
        point_at <0.5,0.5,0.0>
}

// ground -----------------------------------------------------------------
//---------------------------------<<< settings of squared plane dimensions
#declare RasterScale = 0.10;
#declare RasterHalfLine  = 0.0125;
#declare RasterHalfLineZ = 0.0125;
//-------------------------------------------------------------------------
#macro Raster(RScale, HLine)
       pigment{ gradient x scale RScale
                color_map{[0.000   color rgbt<1,1,1,0>*0.8]
                          [0+HLine color rgbt<1,1,1,0>*0.8]
                          [0+HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,1>]
                          [1-HLine color rgbt<1,1,1,0>*0.8]
                          [1.000   color rgbt<1,1,1,0>*0.8]} }
 #end// of Raster(RScale, HLine)-macro
//-------------------------------------------------------------------------

// squared plane XY
plane { <0,0,1>, $floor_z    // plane with layered textures
        texture { pigment{checker color White, color Black}
                scale $light_dy*0.5}
      }
plane { <0,-1,0>, $wall_y    // plane with layered textures
        texture { pigment{color White}
                }
        rotate<0,0,0>
      }
plane { <0,0,-1>, $sky_z    // plane with layered textures
        texture { pigment{color Blue}
                }
        rotate<0,0,0>
      }
//------------------------------------------------ end of squared plane XZ
#macro tank_material()
material{
 texture{
  pigment{ rgbf<.98,.98,.98,0.85>*1}
  finish { ambient 0.0
           diffuse 0.15
           reflection 0.2
           specular 0.6
           roughness 0.005
          // phong 1
          // phong_size 400
           reflection { 0.03, 1.0 fresnel on }
        //   conserve_energy
         }
  } // end of texture

  interior{ ior 1.5
            fade_power 1001
            fade_distance 0.5
            fade_color <0.8,0.8,0.8>
          } // end of interior


} // end of material
#end

#macro isosurface_material()
material{
 texture{
  pigment{ rgbf<.98,.98,.98,0.9>*0.95}
  finish { ambient 0.0
           diffuse 0.15
           specular 0.6
           roughness 0.005
           //phong 1
           //phong_size 400
           reflection { 0.2, 1.0 fresnel on }
           conserve_energy
         }
   } // end of texture

  interior{ ior 1.33
            fade_power 1001
            fade_distance 0.5
            fade_color <0.8,0.8,0.8>
        } // end of interior
} // end of material
#end
""")
            pov.write(povSceneTemplate.substitute(look_at_x=look_at[0],
                                                  look_at_y=look_at[1],
                                                  look_at_z=look_at[2],
                                                  cam_x=cam[0],
                                                  cam_y=cam[1],
                                                  cam_z=cam[2],
                                                  light_x=light[0],
                                                  light_y=light[1],
                                                  light_z=light[2],
                                                  light_dx=light_dx,
                                                  light_dy=light_dy,
                                                  floor_z=floor_z,
                                                  wall_y=wall_y,
                                                  sky_z=sky_z
                                                  )
                      )
            pov.close()
