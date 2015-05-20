from collections import defaultdict, OrderedDict
from itertools import product
from proteus.EGeometry import etriple,ecross,enorm

from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
from numpy.linalg import norm

from . import Comm
from .AuxiliaryVariables import AV_base
from .Profiling import logEvent as log
from proteus.MeshTools import triangleVerticesToNormals, tetrahedronVerticesToNormals, getMeshIntersections


class Isosurface(AV_base):
    """Extract isosurfaces(isocontours)"""
    def __init__(self,isosurfaces, activeTime=None, sampleRate=0, fileName='water.pov'):
        """
        Create a set of isosurfaces that will be extracted and serialized
        
        :param isosurfaces: An iterable of "isosurfaces".  Each isosurface
        is specified by a 2-tuple, with the first element in the tuple is a field
        from which to extract isosurfaces, and the second element is an n-tuple
        of isosurface values to extract.
        :param activeTime: If not None, a 2-tuple of start time and end time for
        which the point gauge is active.
        :param sampleRate: The intervals at which samples should be measured.
        Note that this is a rough lower bound, and that the gauge values could
        be computed less frequently depending on the time integrator.  The default
        value of zero computes the gauge values at every time step.
        :param fileName: The name of the file to serialize results to.

        Example:

        phi0 = (isosurfaces=('phi', (0.0,)),
                activeTime=(0, 2.5),
                sampleRate=0.2,
                fileName='water.pov')

        This creates an Isosurfaces object that will extract phi(x,y,z,t) = 0 at simulation time between = 0 and 2.5 with samples
        taken no more frequently than every 0.2 seconds.  Results will be saved to: water.pov.
        """
        self.isosurfaces = isosurfaces
        for f,values in self.isosurfaces:
            for v in values:
                assert v==0.0,"only implemented for 0 isosurface  in 3D for  now"
        self.activeTime = activeTime
        self.sampleRate = sampleRate
        self.fileName = fileName
    def attachModel(self, model, ar):
        """ Attach this isosurface to the given simulation model.
        """
        self.model = model
        self.fieldNames = model.levelModelList[-1].coefficients.variableNames
        self.elementNodesArray = model.levelModelList[-1].mesh.elementNodesArray
        self.nodeArray = model.levelModelList[-1].mesh.nodeArray
        self.num_owned_elements = model.levelModelList[-1].mesh.nElements_global
        self.u = model.levelModelList[-1].u
        self.timeIntegration = model.levelModelList[-1].timeIntegration
        self.nFrames=0
        self.last_output=None
        return self
    def calculate(self):
        """ Extracts current isosourfaces and updates open output files
        """
        time = self.timeIntegration.tLast
        log("Calculate called at time %g" % time)
        # check that gauge is in its active time region
        if self.activeTime is not None and (self.activeTime[0] > time or self.activeTime[1] < time):
            return

        # check that gauge is ready to be sampled again
        if self.last_output is not None and time < self.last_output + self.sampleRate:
            return
        nodes = []
        elements = []
        normals = []
        normal_indices = []
        assert self.elementNodesArray.shape[1] == 4, "Algorithmn is for tetrahedra but cells have nVertices = d" % (elementNodesArray.shape,)
        for field,values in self.isosurfaces:
            for v in values:
                phi = self.u[self.fieldNames.index(field)].dof
                for eN in range(self.elementNodesArray.shape[0]):
                    plus=[]
                    minus=[]
                    zeros=[]
                    eps = 1.0e-8
                    for  i  in range(4):
                        I = self.elementNodesArray[eN,i]
                        if phi[I] > eps:
                            plus.append(I)
                        elif  phi[I] < eps:
                            minus.append(I)
                        else:
                            zeros.append(I)
                    nZeros = len(zeros)
                    nMinus = len(minus)
                    nPlus = len(plus)
                    assert(nZeros+nMinus+nPlus == 4)
                    nN_start = len(nodes)
                    if nMinus == 2 and nPlus == 2:#4 cut edges
                        for  J in minus:
                            for  I in plus:
                                s = -phi[I]/(phi[J] - phi[I])
                                x = s*(self.nodeArray[J] - self.nodeArray[I])+ self.nodeArray[I]
                                nodes.append(x)
                        elements.append([nN_start+j for j in range(3)])
                        elements.append([nN_start+j+1 for j in range(3)])
                        if etriple(self.nodeArray[plus[0]] - nodes[-4],nodes[-3]-nodes[-4],nodes[-2]- nodes[-4]) < 0.0:
                            elements[-2] = [elements[-2][0], elements[-2][2], elements[-2][1]]
                        if etriple(self.nodeArray[plus[1]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                        normal = ecross(nodes[elements[-2][1]]-nodes[elements[-2][0]],nodes[elements[-2][2]]- nodes[elements[-2][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-2])
                        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-1])
                    elif nPlus == 3 and nMinus == 1:#3 cut edges
                        I = minus[0]
                        for  J in plus:
                            s = -phi[I]/(phi[J] - phi[I])
                            x = s*(self.nodeArray[J] - self.nodeArray[I]) + self.nodeArray[I]
                            nodes.append(x)
                        elements.append([nN_start+j for j in range(3)])
                        if etriple(self.nodeArray[minus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0:
                            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                            assert etriple(self.nodeArray[minus[0]] - nodes[-3],nodes[-1]-nodes[-3],nodes[-2]- nodes[-3]) < 0.0
                        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-1])
                    elif nMinus == 3 and nPlus == 1:#3 cut edges
                        I = plus[0]
                        for  J in minus:
                            s = -phi[I]/(phi[J] - phi[I])
                            x = s*(self.nodeArray[J] - self.nodeArray[I])+ self.nodeArray[I]
                            nodes.append(x)
                        elements.append([nN_start+j for j in range(3)])
                        if etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                            assert etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-1]-nodes[-3],nodes[-2]- nodes[-3]) > 0.0
                        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-1])
                    elif nZeros == 1 and ((nPlus == 2 and nMinus == 1) or
                                          (nPlus == 1 and nMinus == 2)): #2 cut edges, 1 vertex lies in plane
                        nodes.append(self.nodeArray[zeros[0]])
                        for  J in minus:
                            for  I in plus:
                                s = -phi[I]/(phi[J] - phi[I])
                                x = s*(self.nodeArray[J] - self.nodeArray[I])+ self.nodeArray[I]
                                nodes.append(x)
                        elements.append([nN_start+j for j in range(3)])
                        if etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                            assert etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0
                        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-1])
                    elif nZeros == 2 and nPlus == 1 and nMinus == 1:#1 cut edge, 2 vertices lie in plane
                        I = plus[0]
                        J = minus[0]
                        s = -phi[I]/(phi[J] - phi[I])
                        x = s*(self.nodeArray[J] - self.nodeArray[I])+self.nodeArray[I]
                        nodes.append(self.nodeArray[I])
                        nodes.append(self.nodeArray[J])
                        nodes.append(x)
                        elements.append([nN_start+j for j in range(3)])
                        if etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                                elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                                assert etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0
                        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-1])
                    elif nZeros == 3:#3 vertices lie in plane
                        for I in zeros:
                            nodes.append(self.nodeArray[I])
                        elements.append([nN_start+j for  j in range(3)])
                        if nPlus == 1:
                            if etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                                elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                                assert etriple(self.nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0
                        else:
                            if etriple(self.nodeArray[minus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0:
                                elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                                assert etriple(self.nodeArray[minus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0
                        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
                        normal /= enorm(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normals.append(normal)
                        normal_indices.append(elements[-1])
                pov = open(field+"_"+`v`+"_%4.4d.pov" % self.nFrames,"w")
                povScene="""
#include "colors.inc"
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
	location <0.5, -3.0, 1.5>
	sky <0.0, 0.0, 5.0>
	up <0, 0, 1>
	right <1.33, 0, 0>
	angle 45.000000
	look_at <0.5, 0.5, 0.5>
}

light_source {
	<1.5,1.5,1.5>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.5,0.5,0.0>
}
light_source {
	<1.5,-0.5,1.5>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.5,0.5,0.0>
}
light_source {
	<-0.5,-0.5,1.5>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.5,0.5,0.0>
}
light_source {
	<-0.5,1.5,1.5>
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
plane { <0,0,1>, 0.0    // plane with layered textures
        texture { pigment{checker color White, color Black}
                scale 0.04}
      }
plane { <0,0,-1>, 5.0    // plane with layered textures
        texture { pigment{color Blue}
                }
        rotate<0,0,0>
      }
//------------------------------------------------ end of squared plane XZ
 object{
 difference{
 box {
    <0.0-0.01,0.0-0.01,0.0-0.01>,  // Near lower left corner
    <1.0+0.01,1.0+0.01,1.0+0.01>   // Far upper right corner
    }
 box {
    <0.0,0.0,0.0>,  // Near lower left corner
    <1.0,1.0,1.0>   // Far upper right corner
    }}
    
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
  }
  
mesh2 {
        vertex_vectors {
                 """
                pov.write(povScene)
                pov.write("%i,\n" % (len(nodes),))
                for n in nodes:
                    pov.write("<%f,%f,%f>,\n" % tuple(n.tolist()))
                pov.write("""        }
                         normal_vectors {
                                 """)
                pov.write("%i,\n" % (len(normals),))
                for n in normals:
                    pov.write("<%f,%f,%f>,\n" % tuple(n))
                pov.write("""        }
                         face_indices {
                                 """)
                pov.write("%i,\n" % (len(elements),))
                for e in elements:
                    pov.write("<%i,%i,%i>,\n" % tuple(e))
                pov.write("""        }
                         normal_indices {
                                 """)
                pov.write("%i,\n" % (len(normal_indices),))
                for ni in normal_indices:
                    pov.write("<%i,%i,%i>,\n" % tuple(ni))
                pov.write("""        }
                    matrix < 1.000000, 0.000000, 0.000000,
                             0.000000, 1.000000, 0.000000,
                             0.000000, 0.000000, 1.000000,
                             0.000000, 0.000000, 0.000000 >
                        material{
                         texture{
                          pigment{ rgbf<.98,.98,.98,0.9>*0.95}
                          //normal { ripples 1.35 scale 0.0 turbulence 0.3 translate<-0.05,0,0> rotate<0,0,-20>} 
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
                }
                """)
                pov.close()
        self.nFrames+=1
        self.last_output = time
        self.nodes = np.array(nodes)
        self.elements = np.array(elements)
