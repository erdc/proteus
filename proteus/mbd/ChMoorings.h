//#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <string> 
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemDEM.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono_fea/ChElementBeamANCF.h"
#include "chrono_fea/ChElementCableANCF.h"
#include "chrono_fea/ChBeamSection.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChLinkPointPoint.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLoadsBeam.h"
#include "chrono_fea/ChContactSurfaceNodeCloud.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverMINRES.h"


//using namespace std;
using namespace chrono;
using namespace chrono::fea;



class cppCable {
public:
	ChSystemDEM& system;  // global system
	std::shared_ptr<ChMesh> mesh;  // mesh
	int nb_nodes;   // number of nodes along cable
	int nb_elems;   // number of nodes along cable
  double Cd_axial;  // drag coeff in axial direction
  double Cd_normal;  // drag coeff in normal direction
  double Cm_axial;  // added mass coeff in axial direction
  double Cm_normal;  // added mass coeff in normal direction
	std::vector<ChVector<>> mvecs;  // vectors (nodes coordinates)
	std::vector<ChVector<>> mvecs_middle;
	std::vector<ChVector<>> mdirs;  // vectors (nodes coordinates)
	double d, rho, E, length;  // diameter, density, Young's modulus, length of cable
  double A0; // unstretched diameter
	double L0 = 0;  // initial length along cable
  std::string beam_type;
	std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodes;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChElementBeamANCF>> elems;  // array of elements */
	std::vector<std::shared_ptr<ChElementCableANCF>> elemsCableANCF;  // array of elements */
	std::shared_ptr<ChMaterialBeamANCF> mmaterial_cable;  // cable material
	std::vector<std::shared_ptr<ChElementCableANCF>> elems_cable;  // array of elements */
	std::shared_ptr<ChBeamSectionCable> msection_cable;  // cable material
	std::vector<double> elems_length;  // array of elements
	std::vector<ChVector<>> fluid_velocity;
	std::vector<ChVector<>> fluid_acceleration;
	std::vector<double> fluid_density;
  std::vector<double> nodes_density; // density of (cable-fluid) at nodes
	cppCable(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, double length,
           int nb_elems, double d, double rho, double E, double L0, std::string beam_type);  // constructor
	void setFluidVelocityAtNodes(std::vector<ChVector<>> vel);
	void setFluidAccelerationAtNodes(std::vector<ChVector<>> acc);
  void setFluidDensityAtNodes(std::vector<double> vof);
	std::vector<std::shared_ptr<ChVector<double>>> getNodalPositions();
	std::vector<ChVector<>> forces_drag;
	std::vector<ChVector<>> forces_addedmass;
  std::vector<std::shared_ptr<ChLoadBeamWrenchDistributed>> elems_loads_distributed;
  std::vector<std::shared_ptr<ChLoadBeamWrench>> elems_loads;
	void buildVectors();  // builds location vectors for the nodes
	void buildNodes();  // builds the nodes for the mesh
	void buildNodesBeamANCF();  // builds the nodes for the mesh
	void buildNodesCableANCF();  // builds the nodes for the mesh
	void buildMaterials();  // builds the material to use for elements
	void buildMaterialsBeamANCF();  // builds the material to use for elements
	void buildMaterialsCableANCF();  // builds the material to use for elements
	void buildElements();  // builds the elements for the mesh
	void buildElementsBeamANCF();  // builds the elements for the mesh
	void buildElementsCableANCF();  // builds the elements for the mesh
	void buildMesh();  // builds the mesh
	void buildMeshBeamANCF();  // builds the mesh
	void buildMeshCableANCF();  // builds the mesh
	void setDragForce();  // calculates the drag force per nodes
	void setAddedMassForce();  // calculates the added mass force per nodes
  void applyForces();
  void addNodestoContactCloud(std::shared_ptr<ChContactSurfaceNodeCloud> cloud);
  void setDragCoefficients(double axial, double normal);
  void setAddedMassCoefficients(double axial, double normal);
};

class cppMultiSegmentedCable {
public:
	ChSystemDEM& system;  // global system
  std::string beam_type;
	std::shared_ptr<ChMaterialSurfaceDEM> mysurfmaterial;
	std::shared_ptr<ChMesh> mesh;  // mesh
	std::shared_ptr<ChBody> fairleadd;
	std::shared_ptr<ChLinkPointFrame> fairlead2;
	std::vector<int> nb_nodes;   // number of nodes along cable
	std::vector<int> nb_elems;   // number of nodes along cable
	std::vector<ChVector<>> mvecs;  // vectors (nodes coordinates)
	std::vector<std::shared_ptr<cppCable>> cables;
	std::shared_ptr<ChMaterialSurfaceDEM> contact_material;  // mesh
	std::vector<double>  d;
	std::vector<double> rho;
	std::vector<double> E;
	std::vector<double> length;  // diameter, density, Young's modulus, length of cable
	std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodes;  // array nodes coordinates and direction
	std::vector<ChVector<>> fluid_velocity;
	std::vector<ChVector<>> fluid_acceleration;
	std::vector<double> fluid_density;
	std::vector<std::shared_ptr<ChElementBeamANCF>> elems;  // array of elements */
	std::vector<std::shared_ptr<ChElementCableANCF>> elemsCableANCF;  // array of elements */
  std::shared_ptr<ChLinkPointFrame> constraint_front;
  std::shared_ptr<ChLinkPointFrame> constraint_back;
  std::shared_ptr<ChBody> body_back;
  std::shared_ptr<ChBody> body_front;
  bool nodes_built;
	cppMultiSegmentedCable(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, std::vector<double> length,
                         std::vector<int> nb_nodes, std::vector<double> d, std::vector<double> rho, std::vector<double> E, std::string beam_type);
	void setFluidVelocityAtNodes(std::vector<ChVector<>> vel);
	void setFluidAccelerationAtNodes(std::vector<ChVector<>> vel);
    void setFluidDensityAtNodes(std::vector<double> vof);
	void updateDragForces();
	void updateAddedMassForces();
  void applyForces();
	std::vector<std::shared_ptr<ChVector<double>>> getNodalPositions();
  void buildNodes();
	void buildCable();  // builds the multi-segmented cable
	void getForceFairlead();

	void attachBackNodeToBody(std::shared_ptr<ChBody> body);
	void attachFrontNodeToBody(std::shared_ptr<ChBody> body);
  void setContactMaterial(std::shared_ptr<ChMaterialSurfaceDEM> material);
  void buildNodesCloud();
  ChVector<> getTensionElement(int i);
};

class cppMesh {
public:
  ChSystemDEM& system;
	std::shared_ptr<ChMesh> mesh;
	cppMesh(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh);
  void SetAutomaticGravity(bool val);
};

cppMesh::cppMesh(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh) :
  system(system),
	mesh(mesh)
{
  system.Add(mesh);
};

void cppMesh::SetAutomaticGravity(bool val) {
  mesh->SetAutomaticGravity(val);
};


 
cppMultiSegmentedCable::cppMultiSegmentedCable(
	ChSystemDEM& system,
	std::shared_ptr<ChMesh> mesh,
	std::vector<double> length,
	std::vector<int> nb_elems,
	std::vector<double> d,
	std::vector<double> rho,
	std::vector<double> E,
  std::string beam_type = "BeamANCF") :
	system(system),
	mesh(mesh),
	length(length),
	nb_elems(nb_elems),
	d(d),
	rho(rho),
  E(E),
  beam_type(beam_type)
{
  GetLog() << "Initialize cable\n";
  nodes_built = false;
	std::shared_ptr<cppCable> segment;
	double L0 = 0;
	for (int i = 0; i < length.size(); ++i) {
    GetLog() << "Building cable " << i << "\n" ;
		segment = std::make_shared<cppCable>(system, mesh, length[i], nb_elems[i], d[i], rho[i], E[i], L0, beam_type);
		cables.push_back(segment);
		L0 = L0 + length[i];
	}
  GetLog() << "Finished Initializing\n" ;
}

void cppMultiSegmentedCable::buildNodes() {
	nodes.clear();
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->buildNodes();
    nodes.insert(nodes.end(), cables[i]->nodes.begin(), cables[i]->nodes.end());
  }
  nodes_built = true;
}

void cppMultiSegmentedCable::buildCable() {
	/* builds all cable segments and updates their link
	(no duplicate node added to mesh) */
  if (nodes_built == false) {
    nodes.clear();
  }
	elems.clear();
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->buildMaterials();
    if (nodes_built == false) {
      cables[i]->buildNodes();
      nodes.insert(nodes.end(), cables[i]->nodes.begin(), cables[i]->nodes.end());
    }
		cables[i]->buildElements();
    if (beam_type == "BeamANCF") {
      elems.insert(elems.end(), cables[i]->elems.begin(), cables[i]->elems.end());
    }
    else if (beam_type == "CableANCF") {
      elemsCableANCF.insert(elemsCableANCF.end(), cables[i]->elemsCableANCF.begin(), cables[i]->elemsCableANCF.end());
    }
		cables[i]->buildMesh();
		if (i>0) {
			auto con1 = std::make_shared<ChLinkPointPoint>();
			con1->Initialize(cables[i]->nodes.front(), cables[i - 1]->nodes.back());
			system.Add(con1);
		}
	}
  buildNodesCloud();
  fluid_velocity.clear();
  fluid_acceleration.clear();
  fluid_density.clear();
  for (int i = 0; i < nodes.size(); ++i) {
    fluid_velocity.push_back(ChVector<>(0.,0.,0.));
    fluid_acceleration.push_back(ChVector<>(0.,0.,0.));
    fluid_density.push_back(0.);
  }
  setFluidVelocityAtNodes(fluid_velocity);
  setFluidAccelerationAtNodes(fluid_acceleration);
  setFluidDensityAtNodes(fluid_density);
}

void cppMultiSegmentedCable::setFluidAccelerationAtNodes(std::vector<ChVector<>> acc) {
	fluid_acceleration = acc;
	int node_nb = 0;
	for (int i = 0; i < cables.size(); ++i) {
		std::vector<ChVector<>> fluid_acc(fluid_acceleration.begin() + node_nb, fluid_acceleration.begin() + node_nb + cables[i]->nodes.size());
		cables[i]->setFluidAccelerationAtNodes(fluid_acc);
		node_nb += cables[i]->nodes.size();
	}
}

void cppMultiSegmentedCable::setFluidVelocityAtNodes(std::vector<ChVector<>> vel) {
	fluid_velocity = vel;
	int node_nb = 0;
	for (int i = 0; i < cables.size(); ++i) {
		std::vector<ChVector<>> fluid_vel(fluid_velocity.begin() + node_nb, fluid_velocity.begin() + node_nb + cables[i]->nodes.size());
		cables[i]->setFluidVelocityAtNodes(fluid_vel);
		node_nb += cables[i]->nodes.size();
	}
}

void cppMultiSegmentedCable::setFluidDensityAtNodes(std::vector<double> vel) {
	fluid_density = vel;
	int node_nb = 0;
	for (int i = 0; i < cables.size(); ++i) {
		std::vector<double> fluid_dens(fluid_density.begin() + node_nb, fluid_density.begin() + node_nb + cables[i]->nodes.size());
		cables[i]->setFluidDensityAtNodes(fluid_dens);
		node_nb += cables[i]->nodes.size();
	}
}


void cppMultiSegmentedCable::updateDragForces() {
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->setDragForce();
	};
}

void cppMultiSegmentedCable::updateAddedMassForces() {
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->setAddedMassForce();
	};
}

void cppMultiSegmentedCable::applyForces() {
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->applyForces();
	};
}

ChVector<> cppMultiSegmentedCable::getTensionElement(int i) {

  auto mat = ChMatrix<>();
  auto force = ChVector<>();
  auto torque = ChVector<>();
  double eta = 0.6;
  if (beam_type == "BeamANCF") {
    elems[i]->EvaluateSectionForceTorque(eta,
                                         mat,
                                         force,
                                         torque);

  }
  else if (beam_type == "CableANCF") {
    elemsCableANCF[i]->EvaluateSectionForceTorque(eta,
                                         mat,
                                         force,
                                         torque);
  }
  GetLog() << force;
  GetLog() << torque;
  GetLog() << mat;
  return force;
}

std::vector<std::shared_ptr<ChVector<double>>> cppMultiSegmentedCable::getNodalPositions() {
	std::vector<std::shared_ptr<ChVector<double>>> nodal_positions;
	for (int i = 0; i < nodes.size(); ++i) {
		auto pos = nodes[i]->GetPos();
		double x = pos.x();
		double y = pos.y();
		double z = pos.z();
		auto nodal_position = std::make_shared<ChVector<double>>(x, y, z);
		nodal_positions.push_back(nodal_position);
	}
	return nodal_positions;
}

void cppMultiSegmentedCable::attachBackNodeToBody(std::shared_ptr<ChBody> body) {
	auto constraint = std::make_shared<ChLinkPointFrame>();
	constraint->Initialize(nodes.back(), body);
	system.Add(constraint);
  constraint_back = constraint;
  body_back = body;
};

void cppMultiSegmentedCable::attachFrontNodeToBody(std::shared_ptr<ChBody> body) {
	auto constraint = std::make_shared<ChLinkPointFrame>();
	constraint->Initialize(nodes.front(), body);
	system.Add(constraint);
  constraint_front = constraint;
  body_front = body;
};

void cppMultiSegmentedCable::setContactMaterial(std::shared_ptr<ChMaterialSurfaceDEM> material) {
  contact_material = material;
};

void cppMultiSegmentedCable::buildNodesCloud() {
  if (contact_material) {
    auto contact_cloud = std::make_shared<ChContactSurfaceNodeCloud>();
    mesh->AddContactSurface(contact_cloud);
    // Use DEM surface material properties
    contact_cloud->SetMaterialSurface(contact_material);
    // add cable nodes to cloud
    for (int i = 0; i < cables.size(); ++i) {
      cables[i]->addNodestoContactCloud(contact_cloud);
  }
 }
};


cppCable::cppCable(
	ChSystemDEM& system, // system in which the cable belong
	std::shared_ptr<ChMesh> mesh, // mesh of the cable
	double length, // length of cable
	int nb_elems,  // number of nodes along cable
	double d,  // diameter of cable
	double rho,   // density of cable (kg/m3)
	double E, // Young's modulus
	double L0 = 0,
  std::string beam_type = "BeamANCF"
) :
	system(system),
	mesh(mesh),
	length(length),
	nb_elems(nb_elems),
	d(d),
	rho(rho),
	E(E),
    L0(L0),
    beam_type(beam_type)
{
  Cd_axial = 1.15;  // studless chain
  Cd_normal = 1.4;  // studless chain
  Cm_axial = 0.5;  //studless chain
  Cm_normal = 1.;  // studless chain
  A0 = d*d/4*M_PI;
	// TO CHANGE !!!
	//buildMaterials();
	//buildVectors();
	//buildCable(nb_nodes);
	//applyConstraints();
	//double g = { length };
}

void cppCable::buildMaterials() {
  if (beam_type == "BeamANCF") {
    buildMaterialsBeamANCF();
      }
  else if (beam_type == "CableANCF") {
    buildMaterialsCableANCF();
      }
}
void cppCable::buildMaterialsBeamANCF() {
	// make material characteristics of cable
    auto nu = ChVector<>(0.3, 0.3, 0.3);
    double E2 = E / nu.y()*nu.x();
    auto EE = ChVector<>(E, E2, E2);
    double G = EE.z() / (2 * (1 + nu.z()));
    /* double G2 = G; */
    double G2 = 1e-6;
    auto GG = ChVector<>(G, G2, G2);
    double k_rect = 10 * (1 + nu.x()) / (12 + 11 * nu.x()); // rectangular cross-section
    double k_circ = 6 * (1 + nu.x()) / (7 + 6 * nu.x());
    //mmaterial_cable = std::make_shared<ChMaterialBeamANCF>(rho, EE, nu, GG, 0.84, 0.84);
    mmaterial_cable = std::make_shared<ChMaterialBeamANCF>(rho, EE, nu, GG, k_circ, k_circ);
}

void cppCable::buildMaterialsCableANCF() {
    msection_cable = std::make_shared<ChBeamSectionCable>();
    msection_cable->SetDiameter(d);
    msection_cable->SetYoungModulus(E);
    msection_cable->SetDensity(rho);
    msection_cable->SetI(1e-6);
}

void cppCable::buildNodes() {
  GetLog() << "Building nodes\n";
  GetLog() << beam_type << "\n";
  if (beam_type == "BeamANCF") {
    GetLog() << "Building nodes beam\n";
    buildNodesBeamANCF();
  }
  else if (beam_type == "CableANCF") {
    GetLog() << "Building nodes cable\n";
    buildNodesCableANCF();
  }
  GetLog() << "Building nodes finished\n";
}

void cppCable::buildNodesCableANCF() {
	nodes.clear();
  nb_nodes = nb_elems+1;
	ChVector<> dir;  // direction of node
	ChVector<> norm1; // normal of node direction
	ChVector<> norm2; // normal of node direction and norm1
	ChCoordsys<> coordsys;  // coordinate system
	ChVector<> ref = ChVector<>(1., 0., 0.);
	std::shared_ptr<ChNodeFEAxyzDD> node;
	// first node
	dir = (mvecs[1] - mvecs[0]);
	dir.Normalize();
	node = std::make_shared<ChNodeFEAxyzDD>(mvecs[0], dir, ref);
	nodes.push_back(node);
	// other nodes
	for (int i = 1; i < mvecs.size() - 1; ++i) {
		dir = mvecs[i + 1] - mvecs[i - 1];
		dir.Normalize();
		node = std::make_shared<ChNodeFEAxyzDD>(mvecs[i], dir, ref);
		nodes.push_back(node);
	}  // last node
	dir = mvecs[nb_nodes - 1] - mvecs[nb_nodes - 2];
	dir.Normalize();
	node = std::make_shared<ChNodeFEAxyzDD>(mvecs[mvecs.size() - 1], dir, ref);
	nodes.push_back(node);
  GetLog() << "NODES BUILT: " << nodes.size() << "\n";
}

void cppCable::buildNodesBeamANCF() {
	nodes.clear();
  nb_nodes = 2*nb_elems+1;
	ChVector<> dir;  // direction of node
	ChVector<> norm1; // normal of node direction
	ChVector<> norm2; // normal of node direction and norm1
	ChCoordsys<> coordsys;  // coordinate system
	ChVector<> ref = ChVector<>(1., 0., 0.);
	std::shared_ptr<ChNodeFEAxyzDD> node;
	// first node
	dir = (mvecs[1] - mvecs[0]);
	dir.Normalize();
	//plane = dir.x()*(x - mvecs[0].x()) + dir.y()*(y - mvecs[0].y()) + dir.z()*(z - mvecs[0].z());
	if (dir.x() == 1) {
    ref = ChVector<>(0., 0., -1.);
		ref.Normalize();
	}
  else if (dir.x() == -1) {
    ref = ChVector<>(0., 0., 1.);
		ref.Normalize();
  }
	else {
		ref = ChVector<>(1., 0., 0.);
	}
	norm1 = dir % ref;
	norm1.Normalize();
	norm2 = dir % norm1;
	norm2.Normalize();
	//norm2 = dir % norm1;
	node = std::make_shared<ChNodeFEAxyzDD>(mvecs[0], norm1, norm2);
	nodes.push_back(node);
	// other nodes
	for (int i = 1; i < mvecs.size() - 1; ++i) {
		dir = mvecs[i + 1] - mvecs[i - 1];
		dir.Normalize();
    if (dir.x() == 1) {
      ref = ChVector<>(0., 0., -1.);
      ref.Normalize();
    }
    else if (dir.x() == -1) {
      ref = ChVector<>(0., 0., 1.);
      ref.Normalize();
    }
		else {
			ref = ChVector<>(1., 0., 0.);
		}
		norm1 = dir % ref;
		norm1.Normalize();
		norm2 = dir % norm1;
		norm2.Normalize();
		node = std::make_shared<ChNodeFEAxyzDD>(mvecs[i], norm1, norm2);
		nodes.push_back(node);
	}  // last node
	dir = mvecs[nb_nodes - 1] - mvecs[nb_nodes - 2];
	dir.Normalize();
  if (dir.x() == 1) {
    ref = ChVector<>(0., 0., -1.);
    ref.Normalize();
  }
  else if (dir.x() == -1) {
    ref = ChVector<>(0., 0., 1.);
    ref.Normalize();
  }
	else {
		ref = ChVector<>(1., 0., 0.);
	}
	norm1 = dir % ref;
	norm1.Normalize();
	norm2 = dir % norm1;
	norm2.Normalize();
	node = std::make_shared<ChNodeFEAxyzDD>(mvecs[mvecs.size() - 1], norm1, norm2);
	nodes.push_back(node);
}

void cppCable::buildElements() {
  GetLog() << "Building elements";
  if (beam_type == "BeamANCF") {
    GetLog() << "Building Beam elements";
    buildElementsBeamANCF();
      }
  else if (beam_type == "CableANCF") {
    GetLog() << "Building Cable elements";
    buildElementsCableANCF();
      }
  GetLog() << "Finished building elements\n";
}

void cppCable::buildElementsCableANCF() {
	auto loadcontainer = std::make_shared<ChLoadContainer>();
	system.Add(loadcontainer);
	// build elements
	elemsCableANCF.clear();
  elems_loads_distributed.clear();
  elems_loads.clear();
  nb_nodes = nb_elems*2+1;
	for (int i = 0; i < nodes.size() - 1; ++i) {
			auto element = std::make_shared<ChElementCableANCF>();
      auto load_distributed = std::make_shared<ChLoadBeamWrenchDistributed>(element);
      auto load = std::make_shared<ChLoadBeamWrench>(element);
      loadcontainer->Add(load_distributed);
      loadcontainer->Add(load);
			elemsCableANCF.push_back(element);
      elems_loads_distributed.push_back(load_distributed);
      elems_loads.push_back(load);
			element->SetNodes(nodes[i], nodes[i + 1]);
      element->SetSection(msection_cable);
		}
	}

void cppCable::buildElementsBeamANCF() {
	auto loadcontainer = std::make_shared<ChLoadContainer>();
	system.Add(loadcontainer);
	// build elements
	elems.clear();
  elems_loads_distributed.clear();
  elems_loads.clear();
  nb_nodes = nb_elems*2+1;
	for (int i = 1; i < nodes.size() - 1; ++i) {
		if (i % 2 != 0) {
			auto element = std::make_shared<ChElementBeamANCF>();
      auto load_distributed = std::make_shared<ChLoadBeamWrenchDistributed>(element);
      auto load = std::make_shared<ChLoadBeamWrench>(element);
      loadcontainer->Add(load_distributed);
      loadcontainer->Add(load);
			elems.push_back(element);
          elems_loads_distributed.push_back(load_distributed);
          elems_loads.push_back(load);
			element->SetNodes(nodes[i - 1], nodes[i + 1], nodes[i]);
      double elem_length = (nodes[i]->GetPos()-nodes[i-1]->GetPos()).Length()+(nodes[i+1]->GetPos()-nodes[i]->GetPos()).Length();
      double d2 = sqrt(d*d-(d*d-CH_C_PI*d*d/4.));
      d2 = sqrt(d2*d2/2.);  //dividing by 2 because cable thickness seems to be only half of the total thickness in chrono
			element->SetDimensions(elem_length, d2, d2);
			element->SetMaterial(mmaterial_cable);
			element->SetAlphaDamp(0.0004);
			element->SetGravityOn(true);
			element->SetStrainFormulation(ChElementBeamANCF::StrainFormulation::CMNoPoisson);
      element->SetupInitial(&system);
		}
	}
  GetLog() << elems.size() << " Elements built\n";
}

void cppCable::buildMesh() {
  if (beam_type == "BeamANCF") {
    buildMeshBeamANCF();
  }
  else if (beam_type == "CableANCF") {
    buildMeshCableANCF();
  }
}

void cppCable::buildMeshCableANCF() {
	// build the mesh (nodes and elements)
	auto node = elemsCableANCF[0]->GetNodeA();
	mesh->AddNode(node);
	for (int i = 0; i < elemsCableANCF.size(); ++i) {
		auto node = elemsCableANCF[i]->GetNodeB();
		mesh->AddNode(node);
		mesh->AddElement(elemsCableANCF[i]);
	}
}
void cppCable::buildMeshBeamANCF() {
	// build the mesh (nodes and elements)
	auto node = elems[0]->GetNodeA();
	mesh->AddNode(node);
	for (int i = 0; i < elems.size(); ++i) {
		auto node = elems[i]->GetNodeB();
		mesh->AddNode(node);
		node = elems[i]->GetNodeC();
		mesh->AddNode(node);
		mesh->AddElement(elems[i]);
	}
}

void cppCable::setFluidAccelerationAtNodes(std::vector<ChVector<>> acc) {
	fluid_acceleration = acc;
}

void cppCable::setFluidVelocityAtNodes(std::vector<ChVector<>> vel) {
	fluid_velocity = vel;
}

void cppCable::setFluidDensityAtNodes(std::vector<double> dens) {
	fluid_density = dens;
}

void cppCable::setDragCoefficients(double axial, double normal) {
  Cd_axial = axial;
	Cd_normal = normal;
}


void cppCable::setAddedMassCoefficients(double axial, double normal) {
  Cm_axial = axial;
	Cm_normal = normal;
}

void cppCable::setDragForce() {
    /*
     * setFluidVelocityAtNodes and setFluidDensityAtNodes
     * must be called before this function
     */
	ChVector<> u_ch;  // velocity from chrono
	ChVector<> u_prot;  // velocity from proteus
	ChVector<> u_rel;  // relative velocity of node with surrounding fluid
	ChVector<> t_dir;  // tangent at node
	ChVector<> Fd_a;  // axial (tangential) drag force
	ChVector<> Fd_n;  // normal(transversal) drag force
	ChVector<> Fd;  // total drag force
  ChVector<> Va;
  ChVector<> Vn;
  double rho_f;
	// clear current drag forces
	forces_drag.clear();
	double length_elem = length / (nb_nodes - 1);
	for (int i = 0; i < nodes.size(); ++i) {
		u_ch = nodes[i]->GetPos_dt();
		// get velocity u_prot from proteus // TO CHANGE !!
		double ux_prot = fluid_velocity[i][0];
		double uy_prot = fluid_velocity[i][1];
		double uz_prot = fluid_velocity[i][2];
		u_prot = ChVector<>(ux_prot, uy_prot, uz_prot);
		u_rel = u_prot - u_ch;
		// CAREFUL HERE: ChBeamElementANCF, GetD() does not give direction but normal
		t_dir = nodes[i]->GetD() % nodes[i]->GetDD();
		// transverse drag coefficient
		/* double C_dn = 0.;  // transversal drag coeff */
		/* double C_dt = 0.;  // tangential drag coeff */
    rho_f = fluid_density[i];
    Va = u_rel^t_dir*t_dir;
    Vn = u_rel-Va;
    Fd_a = 0.5*rho_f*Cd_axial*d*Va.Length()*Va;//(force per unit length)
    Fd_n = 0.5*rho_f*Cd_normal*M_PI*d*Vn.Length()*Vn;//(force per unit length)
		Fd = Fd_a + Fd_n;
		forces_drag.push_back(Fd);
	}
}


void cppCable::setAddedMassForce() {
    /*
     * setFluidVelocityAtNodes and setFluidDensityAtNodes
     * must be called before this function
     */
	ChVector<> a_ch;  // acceleration from chrono
	ChVector<> a_prot;  // acceleration from proteus
	ChVector<> a_rel;  // relative acceleration of node with surrounding fluid
	ChVector<> t_dir;  // tangent at node
	ChVector<> Fm_a;  // axial (tangential) added mass force
	ChVector<> Fm_n;  // normal(transversal) added mass force
	ChVector<> Fm;  // total added mass force
    ChVector<> Va;
    ChVector<> Vn;
    double rho_f;
	// clear current drag forces
	forces_drag.clear();
	double length_elem = length / (nb_nodes - 1);
	for (int i = 0; i < nodes.size(); ++i) {
		a_ch = nodes[i]->GetPos_dtdt();
		// get velocity u_prot from proteus // TO CHANGE !!
		double ax_prot = fluid_acceleration[i][0];
		double ay_prot = fluid_acceleration[i][1];
		double az_prot = fluid_acceleration[i][2];
		a_prot = ChVector<>(ax_prot, ay_prot, az_prot);
		a_rel = a_prot - a_ch;
		// CAREFUL HERE: ChBeamElementANCF, GetD() does not give direction but normal
		t_dir = nodes[i]->GetD() % nodes[i]->GetDD();
		// transverse drag coefficient
		/* double C_dn = 0.;  // transversal drag coeff */
		/* double C_dt = 0.;  // tangential drag coeff */
    rho_f = fluid_density[i];
    Va = a_rel^t_dir*t_dir;
    Vn = a_rel-Va;
    Fm_a = rho_f*Cm_axial*M_PI*d*d/4.*Va;//(force per unit length)
    Fm_n = rho_f*Cm_normal*M_PI*d*d/4.*Vn;//(force per unit length)
		Fm = Fm_a + Fm_n;
		forces_addedmass.push_back(Fm);
	}
}

void cppCable::applyForces() {
  ChVector<> F_drag;  // drag force per unit length
  ChVector<> F_buoyancy; // buoyancy force per unit length
  ChVector<> F_buoyancy2; // buoyancy force
  ChVector<> F_total;  // total force per unit length
  ChVector<> F_addedmass; // buoyancy force per unit length
  double rho_f;  // density of fluid around cable element
  double mass_f;  // mass of fluid displaced by cable element
  for (int i = 1; i < nodes.size()-1; ++i) {
      if (i % 2 != 0) {
      F_drag = (forces_drag[i-1]+forces_drag[i]+forces_drag[i+1])/3;
      F_addedmass = (forces_addedmass[i-1]+forces_addedmass[i]+forces_addedmass[i+1])/3;
      F_total = F_drag+F_addedmass;//+F_buoyancy;//+F_fluid
      elems_loads_distributed[i/2]->loader.SetForcePerUnit(F_total);
      // buoyancy
      rho_f = -(fluid_density[i-1]+fluid_density[i]+fluid_density[i+1])/3.;
      mass_f = rho_f*A0*elems[i/2]->GetRestLength();
      F_buoyancy = -mass_f*system.Get_G_acc();
      elems_loads[i/2]->loader.SetForce(F_buoyancy);
        }
    }
  };

void cppCable::addNodestoContactCloud(std::shared_ptr<ChContactSurfaceNodeCloud> cloud) {
  for (int i = 0; i < nodes.size(); ++i) {
    cloud->AddNode(nodes[i], d);
      }
};

cppMultiSegmentedCable * newMoorings(
	ChSystemDEM& system,
	std::shared_ptr<ChMesh> mesh,
	std::vector<double> length,
	std::vector<int> nb_elems,
	std::vector<double> d,
	std::vector<double> rho,
	std::vector<double> E,
  std::string beam_type)
{
	return new cppMultiSegmentedCable(system, mesh, length, nb_elems, d, rho, E, beam_type);
}

cppMesh * newMesh(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh) {
	return new cppMesh(system, mesh);
}







class cppSurfaceBoxNodesCloud {
 public:
  ChSystemDEM& system;
  std::shared_ptr<ChMesh> mesh;
  ChVector<> position;
  ChVector<> dimensions;
  std::shared_ptr<ChBodyEasyBox> box;
  std::shared_ptr<ChMaterialSurfaceDEM> material;
  std::shared_ptr<ChContactSurfaceNodeCloud> contact_cloud;
  cppSurfaceBoxNodesCloud(ChSystemDEM& system,
                          std::shared_ptr<ChMesh> mesh,
                          ChVector<> position,
                          ChVector<> dimensions);
  void setNodesSize(double size);
};

cppSurfaceBoxNodesCloud::cppSurfaceBoxNodesCloud(ChSystemDEM& system,
                                                 std::shared_ptr<ChMesh> mesh,
                                                 ChVector<> position,
                                                 ChVector<> dimensions) :
  system(system),
	mesh(mesh),
	position(position),
	dimensions(dimensions)
{
	// Create a surface material to be shared with some objects
	material = std::make_shared<ChMaterialSurfaceDEM>();
	material->SetYoungModulus(2e4);
	material->SetFriction(0.3f);
	material->SetRestitution(0.2f);
	material->SetAdhesion(0);
	contact_cloud = std::make_shared<ChContactSurfaceNodeCloud>();
	mesh->AddContactSurface(contact_cloud);
	// Must use this to 'populate' the contact surface.
	// Use larger point size to match beam section radius
	contact_cloud->AddAllNodes(0.01);

	// Use our DEM surface material properties
	contact_cloud->SetMaterialSurface(material);

	box = std::make_shared<ChBodyEasyBox>(
    dimensions.x(), dimensions.y(), dimensions.z(),  // x,y,z size
		1000,       // density
		true       // collide
		);

	system.Add(box);

	box->SetBodyFixed(true);
	box->SetPos(position);

	// Use our DEM surface material properties
	box->SetMaterialSurface(material);
};

void cppSurfaceBoxNodesCloud::setNodesSize(double size) {
  contact_cloud->AddAllNodes(size);
}

cppSurfaceBoxNodesCloud * newSurfaceBoxNodesCloud(ChSystemDEM& system,
                                                  std::shared_ptr<ChMesh> mesh,
                                                  ChVector<> position,
                                                  ChVector<> dimensions) {
	return new cppSurfaceBoxNodesCloud(system, mesh, position, dimensions);
}
