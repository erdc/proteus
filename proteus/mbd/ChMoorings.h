//#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <string>
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono_fea/ChElementBeamANCF.h"
#include "chrono_fea/ChElementCableANCF.h"
#include "chrono_fea/ChElementBeamEuler.h"
#include "chrono_fea/ChBeamSection.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChLinkPointPoint.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLoadsBeam.h"
#include "chrono_fea/ChContactSurfaceNodeCloud.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/core/ChTransform.h"


//using namespace std;
using namespace chrono;
using namespace chrono::fea;



class MyLoaderTriangular : public ChLoaderUdistributed {
  public:
    // Useful: a constructor that also sets ChLoadable
    ChVector<> Fa;
    ChVector<> Fb;
    MyLoaderTriangular(std::shared_ptr<ChLoadableU> mloadable) : ChLoaderUdistributed(mloadable){
    Fa = ChVector<>(0.,0.,0.);
    Fb = ChVector<>(0.,0.,0.);
    };
    // Compute F=F(u)
    // This is the function that you have to implement. It should return the
    // load at U. For Eulero beams, loads are expected as 6-rows vectors, containing
    // a wrench: forceX, forceY, forceZ, torqueX, torqueY, torqueZ.
    virtual void ComputeF(
        const double U,
        ChVectorDynamic<>& F,
        ChVectorDynamic<>* state_x,
        ChVectorDynamic<>* state_w
        ) {
        double Fy_max = 0.005;
        ChVector<> force = Fa*abs(-1+U)/2.+Fb*(1+U)/2.;
        F.PasteVector(force,0.,0.);  // load, force part; hardwired for brevity
    }
    void SetF(ChVector<> Fa_in, ChVector<> Fb_in) {
    Fa = Fa_in;
    Fb = Fb_in;
    }
    // Needed because inheriting ChLoaderUdistributed. Use 1 because linear load fx.
    virtual int GetIntegrationPointsU() { return 1; }
};
    // Create the load (and handle it with a shared pointer).
    // The ChLoad is a 'container' for your ChLoader.
    // It is created using templates, that is instancing a ChLoad<a_loader_class>()
    //std::shared_ptr<ChLoad<MyLoaderTriangular>> mloadtri(new ChLoad<MyLoaderTriangular>(melementA));
    //mloadcontainer->Add(mloadtri);  // do not forget to add the load to the load container.


class cppCable {
public:
	ChSystemSMC& system;  // global system
	std::shared_ptr<ChMesh> mesh;  // mesh
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
  double Iyy;
  int nb_nodes;
  std::string beam_type;
	std::vector<std::shared_ptr<ChNodeFEAxyzD>> nodes;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodesDD;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChElementCableANCF>> elemsCableANCF;  // array of elements */
	std::vector<std::shared_ptr<ChElementBeamEuler>> elemsBeamEuler;  // array of elements */
	std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodesRot;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChElementCableANCF>> elems_cable;  // array of elements */
	std::shared_ptr<ChBeamSectionCable> msection_cable;  // cable material
	std::shared_ptr<ChBeamSectionAdvanced> msection_advanced;  // cable material
	std::vector<double> elems_length;  // array of elements
	std::vector<ChVector<>> fluid_velocity;
	std::vector<ChVector<>> fluid_acceleration;
	std::vector<double> fluid_density;
  std::vector<double> nodes_density; // density of (cable-fluid) at nodes
	cppCable(ChSystemSMC& system, std::shared_ptr<ChMesh> mesh, double length,
           int nb_elems, double d, double rho, double E, double L0, std::string beam_type);  // constructor
	void setFluidVelocityAtNodes(std::vector<ChVector<>> vel);
	void setFluidAccelerationAtNodes(std::vector<ChVector<>> acc);
  void setFluidDensityAtNodes(std::vector<double> vof);
	std::vector<std::shared_ptr<ChVector<double>>> getNodalPositions();
	std::vector<ChVector<>> forces_drag;
	std::vector<ChVector<>> forces_addedmass;
  std::vector<std::shared_ptr<ChLoadBeamWrenchDistributed>> elems_loads_distributed;
  std::vector<std::shared_ptr<ChLoad<MyLoaderTriangular>>> elems_loads_triangular;
  std::vector<std::shared_ptr<ChLoad<ChLoaderGravity>>> elems_loads_volumetric;
  std::vector<std::shared_ptr<ChLoadBeamWrench>> elems_loads;
	void buildVectors();  // builds location vectors for the nodes
	void buildNodes(bool last_node);  // builds the nodes for the mesh
	void buildNodesBeamEuler(bool last_node);  // builds the nodes for the mesh
	void buildNodesCableANCF(bool last_node);  // builds the nodes for the mesh
	void buildMaterials();  // builds the material to use for elements
	void buildMaterialsCableANCF();  // builds the material to use for elements
  void buildMaterialsBeamEuler();
	void buildElements(bool set_lastnodes);  // builds the elements for the mesh
	void buildElementsBeamEuler(bool set_lastnodes);  // builds the elements for the mesh
	void buildElementsCableANCF(bool set_lastnodes);  // builds the elements for the mesh
	void buildMesh(bool add_lastnode);  // builds the mesh
	void buildMeshBeamEuler(bool add_lastnode);  // builds the mesh
	void buildMeshCableANCF(bool add_lastnode);  // builds the mesh
	void setDragForce();  // calculates the drag force per nodes
	void setAddedMassForce();  // calculates the added mass force per nodes
  void applyForces();
  void addNodestoContactCloud(std::shared_ptr<ChContactSurfaceNodeCloud> cloud);
  void setDragCoefficients(double axial, double normal);
  void setAddedMassCoefficients(double axial, double normal);
  void setIyy(double Iyy_in);
};

class cppMultiSegmentedCable {
public:
	ChSystemSMC& system;  // global system
  std::string beam_type;
	std::shared_ptr<ChMaterialSurfaceSMC> mysurfmaterial;
	std::shared_ptr<ChMesh> mesh;  // mesh
	std::shared_ptr<ChBody> fairleadd;
	std::shared_ptr<ChLinkPointFrame> fairlead2;
	std::vector<int> nb_nodes;   // number of nodes along cable
	std::vector<int> nb_elems;   // number of nodes along cable
	std::vector<ChVector<>> mvecs;  // vectors (nodes coordinates)
	std::vector<std::shared_ptr<cppCable>> cables;
	std::shared_ptr<ChMaterialSurfaceSMC> contact_material;  // mesh
	std::vector<double>  d;
	std::vector<double> rho;
	std::vector<double> E;
	std::vector<double> length;  // diameter, density, Young's modulus, length of cable
	std::vector<std::shared_ptr<ChNodeFEAxyzD>> nodes;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodesDD;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodesRot;  // array nodes coordinates and direction
	std::vector<ChVector<>> fluid_velocity;
	std::vector<ChVector<>> fluid_acceleration;
	std::vector<double> fluid_density;
	std::vector<std::shared_ptr<ChElementCableANCF>> elemsCableANCF;  // array of elements */
	std::vector<std::shared_ptr<ChElementBeamEuler>> elemsBeamEuler;  // array of elements */
  std::shared_ptr<ChLinkPointFrame> constraint_front;
  std::shared_ptr<ChLinkPointFrame> constraint_back;
	std::vector<ChVector<>> forces_drag;
	std::vector<ChVector<>> forces_addedmass;
  std::shared_ptr<ChBody> body_back;
  std::shared_ptr<ChBody> body_front;
  int nb_nodes_tot;
  bool nodes_built;
  bool nodes_chlink;
	cppMultiSegmentedCable(ChSystemSMC& system, std::shared_ptr<ChMesh> mesh, std::vector<double> length,
                         std::vector<int> nb_nodes, std::vector<double> d, std::vector<double> rho, std::vector<double> E, std::string beam_type);
	void setFluidVelocityAtNodes(std::vector<ChVector<>> vel);
	void setFluidAccelerationAtNodes(std::vector<ChVector<>> vel);
  void setFluidDensityAtNodes(std::vector<double> dens);
	void updateDragForces();
	void updateAddedMassForces();
  void applyForces();
	std::vector<std::shared_ptr<ChVector<double>>> getNodalPositions();
  void buildNodes();
	void buildCable();  // builds the multi-segmented cable
	void getForceFairlead();

	void attachBackNodeToBody(std::shared_ptr<ChBody> body);
	void attachFrontNodeToBody(std::shared_ptr<ChBody> body);
  void setContactMaterial(std::shared_ptr<ChMaterialSurfaceSMC> material);
  void buildNodesCloud();
  ChVector<> getTensionElement(int i, double eta);
};

class cppMesh {
public:
  ChSystemSMC& system;
	std::shared_ptr<ChMesh> mesh;
	cppMesh(ChSystemSMC& system, std::shared_ptr<ChMesh> mesh);
  void SetAutomaticGravity(bool val);
};

cppMesh::cppMesh(ChSystemSMC& system, std::shared_ptr<ChMesh> mesh) :
  system(system),
	mesh(mesh)
{
  system.Add(mesh);
};

void cppMesh::SetAutomaticGravity(bool val) {
  mesh->SetAutomaticGravity(val);
};



cppMultiSegmentedCable::cppMultiSegmentedCable(
	ChSystemSMC& system,
	std::shared_ptr<ChMesh> mesh,
	std::vector<double> length,
	std::vector<int> nb_elems,
	std::vector<double> d,
	std::vector<double> rho,
	std::vector<double> E,
  std::string beam_type = "CableANCF") :
	system(system),
	mesh(mesh),
	length(length),
	nb_elems(nb_elems),
	d(d),
	rho(rho),
  E(E),
  beam_type(beam_type)
{
  nodes_built = false;
  nodes_chlink = false;  // build links (true) or link elements directly (false)

	std::shared_ptr<cppCable> segment;
	double L0 = 0;
	for (int i = 0; i < length.size(); ++i) {
		segment = std::make_shared<cppCable>(system, mesh, length[i], nb_elems[i], d[i], rho[i], E[i], L0, beam_type);
		cables.push_back(segment);
		L0 = L0 + length[i];
	}
}

void cppMultiSegmentedCable::buildNodes() {
	nodes.clear();
	nodesRot.clear();
	for (int i = 0; i < cables.size(); ++i) {
    if (beam_type == "BeamEuler") {
        cables[i]->buildNodes(true);
        nodesRot.insert(nodesRot.end(), cables[i]->nodesRot.begin(), cables[i]->nodesRot.end());
        nb_nodes_tot = nodesRot.size();
    }
    else if (beam_type == "CableANCF") {
        if (nodes_chlink == false && i < cables.size()-1) {
            cables[i]->buildNodes(false);
        }
        else {
            cables[i]->buildNodes(true);
        }
        nodes.insert(nodes.end(), cables[i]->nodes.begin(), cables[i]->nodes.end());
        nb_nodes_tot = nodes.size();
        }
    }
  nodes_built = true;
  }

void cppMultiSegmentedCable::buildCable() {
	/* builds all cable segments and updates their link
	(no duplicate node added to mesh) */
  if (nodes_built == false) {
    nodes.clear();
    nodesRot.clear();
    buildNodes();
  }
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->buildMaterials();
    if (i < cables.size()-1 && nodes_chlink == false) {
		cables[i]->buildElements(false);
        if (beam_type == "CableANCF") {
            cables[i]->elemsCableANCF[cables[i]->elemsCableANCF.size()-1]->SetNodes(cables[i]->nodes[cables[i]->nodes.size()-1], cables[i+1]->nodes[0]);
        }
        else if (beam_type == "BeamEuler") {
            cables[i]->elemsBeamEuler[cables[i]->elemsBeamEuler.size()-1]->SetNodes(cables[i]->nodesRot[cables[i]->nodesRot.size()-1], cables[i+1]->nodesRot[0]);
        }
    }
    else {
		cables[i]->buildElements(true);
    }
    if (beam_type == "CableANCF") {
      elemsCableANCF.insert(elemsCableANCF.end(), cables[i]->elemsCableANCF.begin(), cables[i]->elemsCableANCF.end());
    }
    else if (beam_type == "BeamEuler") {
      elemsBeamEuler.insert(elemsBeamEuler.end(), cables[i]->elemsBeamEuler.begin(), cables[i]->elemsBeamEuler.end());
    }
    if (i < cables.size()-1 && nodes_chlink == false) {
        cables[i]->buildMesh(false);
    }
    else {
        cables[i]->buildMesh(true);
    }
    if (nodes_chlink == true) {
        if (i>0) {
            auto con1 = std::make_shared<ChLinkPointPoint>();
            if (beam_type == "BeamEuler") {
                //auto nodeA = cables[i]->nodesRot.front();
                //auto nodeB = cables[i-1]->nodesRot.back();
                //con1->Initialize(nodeA, nodeB);
                //system.Add(con1);
            }
            else if (beam_type == "CableANCF") {
                auto nodeA = cables[i]->nodes.front();
                auto nodeB = cables[i-1]->nodes.back();
                con1->Initialize(nodeA, nodeB);
                system.Add(con1);
            }
        }
    }
	}
  buildNodesCloud();
  fluid_velocity.clear();
  fluid_acceleration.clear();
  fluid_density.clear();
  if (beam_type == "BeamEuler") {
    nb_nodes_tot = nodesRot.size();
  }
  else if (beam_type == "CableANCF") {
    nb_nodes_tot = nodes.size();
  }
  for (int i = 0; i < nb_nodes_tot; ++i) {
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
		std::vector<ChVector<>> fluid_vel(fluid_velocity.begin() + node_nb, fluid_velocity.begin() + node_nb + cables[i]->nb_nodes);
		cables[i]->setFluidVelocityAtNodes(fluid_vel);
		node_nb += cables[i]->nb_nodes;
	}
}

void cppMultiSegmentedCable::setFluidDensityAtNodes(std::vector<double> dens) {
	fluid_density = dens;
	int node_nb = 0;
	for (int i = 0; i < cables.size(); ++i) {
		std::vector<double> fluid_dens(fluid_density.begin() + node_nb, fluid_density.begin() + node_nb + cables[i]->nb_nodes);
		cables[i]->setFluidDensityAtNodes(fluid_dens);
		node_nb += cables[i]->nb_nodes;
	}
}


void cppMultiSegmentedCable::updateDragForces() {
  forces_drag.clear();
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->setDragForce();
    forces_drag.insert(forces_drag.end(), cables[i]->forces_drag.begin(), cables[i]->forces_drag.end());
	};
}

void cppMultiSegmentedCable::updateAddedMassForces() {
  forces_addedmass.clear();
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->setAddedMassForce();
    forces_addedmass.insert(forces_addedmass.end(), cables[i]->forces_addedmass.begin(), cables[i]->forces_addedmass.end());
	};
}

void cppMultiSegmentedCable::applyForces() {
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->applyForces();
	};
}

ChVector<> cppMultiSegmentedCable::getTensionElement(int i, double eta=0.) {

  auto mat = ChMatrix<>();
  auto force = ChVector<>();
  auto torque = ChVector<>();
  if (beam_type == "CableANCF") {
    elemsCableANCF[i]->EvaluateSectionForceTorque(eta,
                                                  mat,
                                                  force,
                                                  torque);
    elemsCableANCF[i]->EvaluateSectionStrain(eta, mat, force);
  }
  else if (beam_type == "BeamEuler") {
    elemsBeamEuler[i]->EvaluateSectionStrain(eta, mat, force);
    /* elemsBeamEuler[i]->EvaluateSectionForceTorque(eta, */
    /*                                               mat, */
    /*                                              force, */
    /*                                               torque); */
  }
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
  if (beam_type == "BeamEuler") {
    //constraint->Initialize(nodesRot.back(), body);
    /* system.Add(constraint); */
    /* constraint_back = constraint; */
    /* body_back = body; */
  }
  else {
    constraint->Initialize(nodes.back(), body);
    system.Add(constraint);
    constraint_back = constraint;
    body_back = body;
  }
};

void cppMultiSegmentedCable::attachFrontNodeToBody(std::shared_ptr<ChBody> body) {
	auto constraint = std::make_shared<ChLinkPointFrame>();
  if (beam_type == "BeamEuler") {
    //constraint->Initialize(nodesRot.front(), body);
    /* system.Add(constraint); */
    /* constraint_front = constraint; */
    /* body_front = body; */
  }
  else if (beam_type == "CableANCF") {
    constraint->Initialize(nodes.front(), body);
    system.Add(constraint);
    constraint_front = constraint;
    body_front = body;
  }
};

void cppMultiSegmentedCable::setContactMaterial(std::shared_ptr<ChMaterialSurfaceSMC> material) {
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
	ChSystemSMC& system, // system in which the cable belong
	std::shared_ptr<ChMesh> mesh, // mesh of the cable
	double length, // length of cable
	int nb_elems,  // number of nodes along cable
	double d,  // diameter of cable
	double rho,   // density of cable (kg/m3)
	double E, // Young's modulus
	double L0 = 0,
  std::string beam_type = "CableANCF"
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
  Iyy = 1e-12;
	// TO CHANGE !!!
	//buildMaterials();
	//buildVectors();
	//buildCable(nb_nodes);
	//applyConstraints();
	//double g = { length };
}

void cppCable::buildMaterials() {
  if (beam_type == "CableANCF") {
    buildMaterialsCableANCF();
      }
  else if (beam_type == "BeamEuler") {
    buildMaterialsBeamEuler();
  }
}
void cppCable::buildMaterialsCableANCF() {
    msection_cable = std::make_shared<ChBeamSectionCable>();
    msection_cable->SetDiameter(d);
    msection_cable->SetYoungModulus(E);
    msection_cable->SetDensity(rho);
    msection_cable->SetI(Iyy);
}

void cppCable::buildMaterialsBeamEuler() {
  msection_advanced = std::make_shared<ChBeamSectionAdvanced>();
  msection_advanced->SetYoungModulus(E);
  msection_advanced->SetGshearModulus(1e-6);
  msection_advanced->SetDensity(rho);
  msection_advanced->SetAsCircularSection(d);
  msection_advanced->SetIyy(Iyy);
  msection_advanced->SetIzz(Iyy);
}

void cppCable::setIyy(double Iyy_in) {
  Iyy = Iyy_in;
}

void cppCable::buildNodes(bool last_node=true) {
  if (beam_type == "CableANCF") {
    buildNodesCableANCF(last_node);
  }
  else if (beam_type == "BeamEuler") {
    buildNodesBeamEuler(last_node);
  }
}

void cppCable::buildNodesBeamEuler(bool last_node) {
	nodesRot.clear();
	ChVector<> dir;  // direction of node
	/* ChQuaternion<> quat;  // direction of node */
	std::shared_ptr<ChNodeFEAxyzrot> node;
	// first node
	dir = (mvecs[1] - mvecs[0]);
	dir.Normalize();
  /* quat = Q_from_AngAxis(dir, 0.); */
	node = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(mvecs[0]));
	nodesRot.push_back(node);
	// other nodes
	for (int i = 1; i < mvecs.size() - 1; ++i) {
		dir = mvecs[i + 1] - mvecs[i - 1];
		dir.Normalize();
    /* quat = Q_from_AngAxis(dir, 0.); */
		node = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>( mvecs[i] ));
		nodesRot.push_back(node);
	}  // last node
    if (last_node == true) {
        nb_nodes = nb_elems+1;
        dir = mvecs[nb_nodes - 1] - mvecs[nb_nodes - 2];
        dir.Normalize();
    /* quat = Q_from_AngAxis(dir, 0.); */
        node = std::make_shared<ChNodeFEAxyzrot>(ChFrame<>(mvecs[mvecs.size() - 1]));
        nodesRot.push_back(node);
    }
    else {
        nb_nodes = nb_elems;
    }
}

void cppCable::buildNodesCableANCF(bool last_node) {
	nodes.clear();
	ChVector<> dir;  // direction of node
	ChVector<> norm1; // normal of node direction
	ChVector<> norm2; // normal of node direction and norm1
	ChCoordsys<> coordsys;  // coordinate system
	ChVector<> ref = ChVector<>(1., 0., 0.);
	std::shared_ptr<ChNodeFEAxyzD> node;
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
    if (last_node == true) {\
        nb_nodes = nb_elems+1;
        dir = mvecs[nb_nodes - 1] - mvecs[nb_nodes - 2];
        dir.Normalize();
        node = std::make_shared<ChNodeFEAxyzDD>(mvecs[mvecs.size() - 1], dir, ref);
        nodes.push_back(node);
    }
    else {
        nb_nodes = nb_elems;
    }
}

void cppCable::buildElements(bool set_lastnodes=true) {
  if (beam_type == "CableANCF") {
    buildElementsCableANCF(set_lastnodes);
      }
  else if (beam_type == "BeamEuler") {
    buildElementsBeamEuler(set_lastnodes);
  }
}

void cppCable::buildElementsCableANCF(bool set_lastnodes) {
	auto loadcontainer = std::make_shared<ChLoadContainer>();
	system.Add(loadcontainer);
	// build elements
	elemsCableANCF.clear();
      elems_loads_distributed.clear();
      elems_loads_volumetric.clear();
      elems_loads_triangular.clear();
      elems_loads.clear();
	for (int i = 0; i < nb_elems; ++i) {
			auto element = std::make_shared<ChElementCableANCF>();
      auto load_distributed = std::make_shared<ChLoadBeamWrenchDistributed>(element);
      auto load = std::make_shared<ChLoadBeamWrench>(element);

      std::shared_ptr<ChLoad<MyLoaderTriangular>> loadtri(new ChLoad<MyLoaderTriangular>(element));
      auto load_volumetric = std::make_shared<ChLoad<ChLoaderGravity>>(element);
      loadcontainer->Add(load_distributed);
      loadcontainer->Add(load);
      loadcontainer->Add(loadtri);  // do not forget to add the load to the load container.
      loadcontainer->Add(load_volumetric);
      elemsCableANCF.push_back(element);
      elems_loads_distributed.push_back(load_distributed);
      elems_loads.push_back(load);
      elems_loads_triangular.push_back(loadtri);
      elems_loads_volumetric.push_back(load_volumetric);
      element->SetSection(msection_cable);
      if (i < nb_elems-1) {
        element->SetNodes(nodes[i], nodes[i + 1]);
      }
      else {
        if (set_lastnodes == true) {
            element->SetNodes(nodes[i], nodes[i + 1]);
        }
      }
    }
	}

void cppCable::buildElementsBeamEuler(bool set_lastnodes) {
	auto loadcontainer = std::make_shared<ChLoadContainer>();
	system.Add(loadcontainer);
	// build elements
	elemsBeamEuler.clear();
  elems_loads_distributed.clear();
  elems_loads.clear();
	for (int i = 0; i < nodesRot.size() - 1; ++i) {
    auto element = std::make_shared<ChElementBeamEuler>();
    auto load_distributed = std::make_shared<ChLoadBeamWrenchDistributed>(element);
    auto load = std::make_shared<ChLoadBeamWrench>(element);
    loadcontainer->Add(load_distributed);
    loadcontainer->Add(load);
    elemsBeamEuler.push_back(element);
    elems_loads_distributed.push_back(load_distributed);
    elems_loads.push_back(load);
    element->SetSection(msection_advanced);
    if (i < nb_elems-1) {
        element->SetNodes(nodesRot[i], nodesRot[i + 1]);
    }
    else {
        if (set_lastnodes == true) {
            element->SetNodes(nodesRot[i], nodesRot[i + 1]);
        }
    }
  }
}

void cppCable::buildMesh(bool add_lastnode=true) {
  if (beam_type == "CableANCF") {
    buildMeshCableANCF(add_lastnode);
  }
  else if (beam_type == "BeamEuler") {
    buildMeshBeamEuler(add_lastnode);
  }
}

void cppCable::buildMeshBeamEuler(bool add_lastnode) {
	// build the mesh (nodes and elements)
	auto node = elemsBeamEuler[0]->GetNodeA();
	mesh->AddNode(node);
	for (int i = 0; i < elemsBeamEuler.size()-1; ++i) {
		auto node = elemsBeamEuler[i]->GetNodeB();
		mesh->AddNode(node);
		mesh->AddElement(elemsBeamEuler[i]);
	}
    if (add_lastnode == true) {
		auto node = elemsBeamEuler[elemsBeamEuler.size()-1]->GetNodeB();
		mesh->AddNode(node);
    }
    mesh->AddElement(elemsBeamEuler[elemsBeamEuler.size()-1]);
}

void cppCable::buildMeshCableANCF(bool add_lastnode) {
	// build the mesh (nodes and elements)
	auto node = elemsCableANCF[0]->GetNodeA();
	mesh->AddNode(node);
	for (int i = 0; i < elemsCableANCF.size()-1; ++i) {
		auto node = elemsCableANCF[i]->GetNodeB();
		mesh->AddElement(elemsCableANCF[i]);
		mesh->AddNode(node);
	}
    if (add_lastnode == true) {
		auto node = elemsCableANCF[elemsCableANCF.size()-1]->GetNodeB();
		mesh->AddNode(node);
    }
    mesh->AddElement(elemsCableANCF[elemsCableANCF.size()-1]);
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
	for (int i = 0; i < nb_nodes; ++i) {
		u_ch = nodes[i]->GetPos_dt();
		// get velocity u_prot from proteus // TO CHANGE !!
		double ux_prot = fluid_velocity[i][0];
		double uy_prot = fluid_velocity[i][1];
		double uz_prot = fluid_velocity[i][2];
		u_prot = ChVector<>(ux_prot, uy_prot, uz_prot);
		u_rel = u_prot - u_ch;
		// CAREFUL HERE: ChBeamElementANCF, GetD() does not give direction but normal
        if (beam_type == "CableANCF") {
            t_dir = nodes[i]->GetD();
        }
        else if (beam_type == "BeamEuler") {
            //t_dir = nodesRot[i]->GetD();
        }
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
	forces_addedmass.clear();
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
        if (beam_type == "CableANCF") {
            t_dir = nodes[i]->GetD();
        }
        else if (beam_type == "BeamEuler") {
            //t_dir = nodesRot[i]->GetD();
        }
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
  ChVector<> F_total;  // total force per unit length
  ChVector<> F_addedmass; // buoyancy force per unit length
  double rho_f;  // density of fluid around cable element
  double lengths = 0;
  for (int i = 0; i < nodes.size()-1; ++i) {
      //F_drag = (forces_drag[i]+forces_drag[i+1])/2;
      //F_addedmass = (forces_addedmass[i]+forces_addedmass[i+1])/2;
      //F_total = F_drag+F_addedmass;
      //elems_loads_distributed[i]->loader.SetForcePerUnit(F_total);
      ChVector<> Fa = forces_drag[i]+forces_addedmass[i];
      ChVector<> Fb = forces_drag[i+1]+forces_addedmass[i+1];
      elems_loads_triangular[i]->loader.SetF(Fa, Fb);
      // buoyancy
      elems_loads_volumetric[i]->loader.Set_G_acc(-fluid_density[i]/rho*system.Get_G_acc()); // remove gravity because it is there on top of it
    }
  };

void cppCable::addNodestoContactCloud(std::shared_ptr<ChContactSurfaceNodeCloud> cloud) {
  for (int i = 0; i < nb_nodes; ++i) {
    if (beam_type == "BeamEuler") {
      cloud->AddNode(nodesRot[i], d);
    }
    else if (beam_type == "CableANCF") {
      cloud->AddNode(nodes[i], d);
    }
  }
};

cppMultiSegmentedCable * newMoorings(
	ChSystemSMC& system,
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

cppMesh * newMesh(ChSystemSMC& system, std::shared_ptr<ChMesh> mesh) {
	return new cppMesh(system, mesh);
}







class cppSurfaceBoxNodesCloud {
 public:
  ChSystemSMC& system;
  std::shared_ptr<ChMesh> mesh;
  ChVector<> position;
  ChVector<> dimensions;
  std::shared_ptr<ChBodyEasyBox> box;
  std::shared_ptr<ChMaterialSurfaceSMC> material;
  std::shared_ptr<ChContactSurfaceNodeCloud> contact_cloud;
  cppSurfaceBoxNodesCloud(ChSystemSMC& system,
                          std::shared_ptr<ChMesh> mesh,
                          ChVector<> position,
                          ChVector<> dimensions);
  void setNodesSize(double size);
};

cppSurfaceBoxNodesCloud::cppSurfaceBoxNodesCloud(ChSystemSMC& system,
                                                 std::shared_ptr<ChMesh> mesh,
                                                 ChVector<> position,
                                                 ChVector<> dimensions) :
  system(system),
	mesh(mesh),
	position(position),
	dimensions(dimensions)
{
	// Create a surface material to be shared with some objects
	material = std::make_shared<ChMaterialSurfaceSMC>();
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

cppSurfaceBoxNodesCloud * newSurfaceBoxNodesCloud(ChSystemSMC& system,
                                                  std::shared_ptr<ChMesh> mesh,
                                                  ChVector<> position,
                                                  ChVector<> dimensions) {
	return new cppSurfaceBoxNodesCloud(system, mesh, position, dimensions);
}

/* class MyLoaderTriangular : public ChLoaderUdistributed { */
/*  public: */
/*   // Useful: a constructor that also sets ChLoadable */
/*  MyLoaderTriangular(std::shared_ptr<ChLoadableU> mloadable) : ChLoaderUdistributed(mloadable){}; */
/*   // Compute F=F(u) */
/*   // This is the function that you have to implement. It should return the */
/*   // load at U. For Eulero beams, loads are expected as 6-rows vectors, containing */
/*   // a wrench: forceX, forceY, forceZ, torqueX, torqueY, torqueZ. */
/*   virtual void ComputeF( */
/*                         const double U,               */
/*                         ChVectorDynamic<>& F,         */
/*                         ChVectorDynamic<>* state_x,   */
/*                         ChVectorDynamic<>* state_w    */
/*                         ) { */
/*     double Fy_max = 0.005; */
/*     ChVector<> force = -(U-1)/2*F1 + (U+1)/2*F2 */
/*     F.PasteVector(ChVector<>) */
/*     F.PasteVector(ChVector<>(0, (((1 + U) / 2) * Fy_max), 0), 0, 0);  // load, force part; hardwired for brevity */
/*     F.PasteVector(ChVector<>(0, 0, 0), 3, 0);  // load, torque part; hardwired for brevity */
/*   } */
