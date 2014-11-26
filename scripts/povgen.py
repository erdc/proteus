import tables
import numpy as np
from proteus.EGeometry import etriple,ecross,enorm
#todo
#set up camera, lighting, ground, and sky from bounding box
#use PLC domain to generate glass container
#add proper script interface
#add support for partitioned domains and multiple timesteps
h5 = tables.openFile("proteus_dem0.h5","r")
pov = open("water.pov","w")

nodeArray = h5.getNode("/nodesSpatial_Domain0")[:]
elementNodesArray = h5.getNode("/elementsSpatial_Domain0")[:]
phi = h5.getNode("/phi10")[:]
ll = (nodeArray[:,0].min(),nodeArray[:,1].min(),nodeArray[:,2].min())
ur = (nodeArray[:,0].max(),nodeArray[:,1].max(),nodeArray[:,2].max())
L = tuple(np.array(ur) - np.array(ll))
center = (nodeArray[:,0].mean(),nodeArray[:,1].mean(),nodeArray[:,2].mean())
print ll,ur,center,L
nodes = []
elements = []
normals = []
normal_indices = []
assert elementNodesArray.shape[1] == 4, "Algorithmn is for tetrahedra but cells have nVertices = d" % (elementNodesArray.shape,)
for eN in range(elementNodesArray.shape[0]):
    plus=[]
    minus=[]
    zeros=[]
    eps = 1.0e-8
    for  i  in range(4):
        I = elementNodesArray[eN,i]
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
                x = s*(nodeArray[J] - nodeArray[I])+ nodeArray[I]
                nodes.append(x)
        elements.append([nN_start+j for j in range(3)])
        elements.append([nN_start+j+1 for j in range(3)])
        if etriple(nodeArray[plus[0]] - nodes[-4],nodes[-3]-nodes[-4],nodes[-2]- nodes[-4]) < 0.0:
            elements[-2] = [elements[-2][0], elements[-2][2], elements[-2][1]]
        if etriple(nodeArray[plus[1]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
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
            x = s*(nodeArray[J] - nodeArray[I]) + nodeArray[I]
            nodes.append(x)
        elements.append([nN_start+j for j in range(3)])
        if etriple(nodeArray[minus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0:
            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
            assert etriple(nodeArray[minus[0]] - nodes[-3],nodes[-1]-nodes[-3],nodes[-2]- nodes[-3]) < 0.0
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
            x = s*(nodeArray[J] - nodeArray[I])+ nodeArray[I]
            nodes.append(x)
        elements.append([nN_start+j for j in range(3)])
        if etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
            assert etriple(nodeArray[plus[0]] - nodes[-3],nodes[-1]-nodes[-3],nodes[-2]- nodes[-3]) > 0.0
        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
        normal /= enorm(normal)
        normals.append(normal)
        normals.append(normal)
        normals.append(normal)
        normal_indices.append(elements[-1])
    elif nZeros == 1 and ((nPlus == 2 and nMinus == 1) or
                          (nPlus == 1 and nMinus == 2)): #2 cut edges, 1 vertex lies in plane
        nodes.append(nodeArray[zeros[0]])
        for  J in minus:
            for  I in plus:
                s = -phi[I]/(phi[J] - phi[I])
                x = s*(nodeArray[J] - nodeArray[I])+ nodeArray[I]
                nodes.append(x)
        elements.append([nN_start+j for j in range(3)])
        if etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
            elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
            assert etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0
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
        x = s*(nodeArray[J] - nodeArray[I])+nodeArray[I]
        nodes.append(nodeArray[I])
        nodes.append(nodeArray[J])
        nodes.append(x)
        elements.append([nN_start+j for j in range(3)])
        if etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                assert etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0
        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
        normal /= enorm(normal)
        normals.append(normal)
        normals.append(normal)
        normals.append(normal)
        normal_indices.append(elements[-1])
    elif nZeros == 3:#3 vertices lie in plane
        for I in zeros:
            nodes.append(nodeArray[I])
        elements.append([nN_start+j for  j in range(3)])
        if nPlus == 1:
            if etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0:
                elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                assert etriple(nodeArray[plus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0
        else:
            if etriple(nodeArray[minus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) > 0.0:
                elements[-1] = [elements[-1][0], elements[-1][2], elements[-1][1]]
                assert etriple(nodeArray[minus[0]] - nodes[-3],nodes[-2]-nodes[-3],nodes[-1]- nodes[-3]) < 0.0
        normal = ecross(nodes[elements[-1][1]]-nodes[elements[-1][0]],nodes[elements[-1][2]]- nodes[elements[-1][0]])
        normal /= enorm(normal)
        normals.append(normal)
        normals.append(normal)
        normals.append(normal)
        normal_indices.append(elements[-1])

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
	location <0.0, -0.75, 0.2>
	sky <0.0, 0.0, 5.0>
	up <0, 0, 1>
	right <1.33, 0, 0>
	angle 45.000000
	look_at <0.0, 0.0, -0.1>
}

light_source {
	<0.2,0.2,0.4>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.0,0.0,0.0>
}
light_source {
	<0.2,-0.2,0.4>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.0,0.0,0.0>
}
light_source {
	<-0.2,-0.2,0.4>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.0,0.0,0.0>
}
light_source {
	<-0.2,0.2,0.4>
	color <0.99980, 0.999800, 0.999800>*2.250000
	spotlight
	point_at <0.0,0.0,0.0>
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
plane { <0,0,1>, -0.11    // plane with layered textures
        texture { pigment{checker color White, color Black}
                scale 0.04}
      }
plane { <0,0,-1>, 1.0    // plane with layered textures
        texture { pigment{color Blue}
                }
        rotate<0,0,0>
      }
//------------------------------------------------ end of squared plane XZ
 object{
 difference{
 box {
    <-0.2-0.01,-0.2-0.01,-0.1>,  // Near lower left corner
    <0.2+0.01,0.2+0.01,0.1>   // Far upper right corner
    }
 box {
    <-0.2,-0.2,-0.1+0.01>,  // Near lower left corner
    <0.2,0.2,0.1+0.01>   // Far upper right corner
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
//material
//{texture{
//pigment { color Blue }
//T_Stone25 scale 2
//    pigment{ rgbf<.93,.95,.98,0.825>*0.99}
//    finish { ambient 0.0 diffuse 0.15
//             reflection{0.1,0.1}
//             specular 0.6 roughness 0.005
//             conserve_energy
//           } // end finish
//  } // end of texture

//    interior{ ior 1.33
//             fade_power 1001
//             fade_distance 0.5
//             fade_color <0.8,0.8,0.8>
//             caustics 0.16
//   } // end of interior
//}
}
""")
h5.close()
pov.close()
