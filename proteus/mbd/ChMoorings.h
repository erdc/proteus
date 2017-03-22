//#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <string> 
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemDEM.h"
#include "chrono/physics/ChLoadContainer.h"
/* #include "chrono/physics/ChBodyEasy.h" */
#include "chrono_fea/ChElementBeamANCF.h"
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
	std::vector<ChVector<>> mvecs;  // vectors (nodes coordinates)
	std::vector<ChVector<>> mvecs_middle;
	std::vector<ChVector<>> mdirs;  // vectors (nodes coordinates)
	std::shared_ptr<ChMaterialBeamANCF> mmaterial_cable;  // cable material
	double d, rho, E, length;  // diameter, density, Young's modulus, length of cable
	double L0 = 0;  // initial length along cable
	std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodes;  // array nodes coordinates and direction
	std::vector<std::shared_ptr<ChElementBeamANCF>> elems;  // array of elements */
	std::vector<double> elems_length;  // array of elements
	cppCable(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, double length,
		int nb_elems, double d, double rho, double E, double L0);  // constructor
	void setFluidVelocityAtNodes(std::vector<ChVector<>> vel);
	std::vector<ChVector<>> fluid_velocity;
	std::vector<std::shared_ptr<ChVector<double>>> getNodalPositions();
	std::vector < std::shared_ptr<ChVector<>>> forces_drag;
	void buildVectors();  // builds location vectors for the nodes
	void buildNodes();  // builds the nodes for the mesh
	void buildMaterials();  // builds the material to use for elements
	void buildElements();  // builds the elements for the mesh
	void buildMesh();  // builds the mesh
	void setDragForce();  // calculates the drag force per nodes
	void getBuoyancyForce();  // calculates buoyancy force
};

class cppMultiSegmentedCable {
public:
	ChSystemDEM& system;  // global system
	std::shared_ptr<ChMaterialSurfaceDEM> mysurfmaterial;
	std::shared_ptr<ChMesh> mesh;  // mesh
	std::shared_ptr<ChBody> fairleadd;
	std::shared_ptr<ChLinkPointFrame> fairlead2;
	std::vector<int> nb_nodes;   // number of nodes along cable
	std::vector<int> nb_elems;   // number of nodes along cable
	std::vector<ChVector<>> mvecs;  // vectors (nodes coordinates)
	std::vector<std::shared_ptr<cppCable>> cables;
	std::vector<double>  d;
	std::vector<double> rho;
	std::vector<double> E;
	std::vector<double> length;  // diameter, density, Young's modulus, length of cable
	std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodes;  // array nodes coordinates and direction
	std::vector<ChVector<>> fluid_velocity;
	std::vector<std::shared_ptr<ChElementBeamANCF>> elems;  // array of elements */
	cppMultiSegmentedCable(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, std::vector<double> length,
		std::vector<int> nb_nodes, std::vector<double> d, std::vector<double> rho, std::vector<double> E);
	void setFluidVelocityAtNodes(std::vector<ChVector<>> vel);
	void updateDragForces();
	std::vector<std::shared_ptr<ChVector<double>>> getNodalPositions();
	void buildCable();  // builds the multi-segmented cable
	void getForceFairlead();

	void attachBackNodeToBody(std::shared_ptr<ChBody> body);
	void attachFrontNodeToBody(std::shared_ptr<ChBody> body);
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
	std::vector<double> E) :
	system(system),
	mesh(mesh),
	length(length),
	nb_elems(nb_elems),
	d(d),
	rho(rho),
	E(E)
{
	std::shared_ptr<cppCable> segment;
	double L0 = 0;
	for (int i = 0; i < length.size(); ++i) {
		segment = std::make_shared<cppCable>(system, mesh, length[i], nb_elems[i], d[i], rho[i], E[i], L0);
		cables.push_back(segment);
		L0 = L0 + length[i];
	}
}

void cppMultiSegmentedCable::buildCable() {
	/* builds all cable segments and updates their link
	(no duplicate node added to mesh) */
	nodes.clear();
	elems.clear();
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->buildMaterials();
		cables[i]->buildNodes();
		nodes.insert(nodes.end(), cables[i]->nodes.begin(), cables[i]->nodes.end());
		cables[i]->buildElements();
		cables[i]->buildMesh();
		if (i>0) {
			auto con1 = std::make_shared<ChLinkPointPoint>();
			con1->Initialize(cables[i]->nodes.front(), cables[i - 1]->nodes.back());
			system.Add(con1);
		}
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

void cppMultiSegmentedCable::updateDragForces() {
	for (int i = 0; i < cables.size(); ++i) {
		cables[i]->setDragForce();
	};
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
};

void cppMultiSegmentedCable::attachFrontNodeToBody(std::shared_ptr<ChBody> body) {
	auto constraint = std::make_shared<ChLinkPointFrame>();
	constraint->Initialize(nodes.front(), body);
	system.Add(constraint);
};

cppCable::cppCable(
	ChSystemDEM& system, // system in which the cable belong
	std::shared_ptr<ChMesh> mesh, // mesh of the cable
	double length, // length of cable
	int nb_elems,  // number of nodes along cable
	double d,  // diameter of cable
	double rho,   // density of cable (kg/m3)
	double E, // Young's modulus
	double L0 = 0
) :
	system(system),
	mesh(mesh),
	length(length),
	nb_elems(nb_elems),
	d(d),
	rho(rho),
	E(E),
	L0(L0)
{
	// TO CHANGE !!!
	//buildMaterials();
	//buildVectors();
	//buildCable(nb_nodes);
	//applyConstraints();
	//double g = { length };
}

void cppCable::buildMaterials() {
	// make material characteristics of cable
	auto nu = ChVector<>(0.3, 0.3, 0.3);
	double E2 = E / nu.y()*nu.x();
	auto EE = ChVector<>(E, E2, E2);
	double G2 = 1e-6;
	auto GG = ChVector<>(EE.z() / (2 * (1 + nu.z())), G2, G2);
	double k_rect = 10 * (1 + nu.x()) / (12 + 11 * nu.x()); // rectangular cross-section
	double k_circ = 6 * (1 + nu.x()) / (7 + 6 * nu.x());
	//mmaterial_cable = std::make_shared<ChMaterialBeamANCF>(rho, EE, nu, GG, 0.84, 0.84);
	mmaterial_cable = std::make_shared<ChMaterialBeamANCF>(rho, EE, nu, GG, k_circ, k_circ);
}

void cppCable::buildVectors() {
	// find coordinates of nodes
	// s = a.sinh(x/a); x = a.arcsinh(s/a); y = a.cosh(x/a)
	double s = 26;
	double a = 6.5256159;
	double x0 = 0;
	// fraction of length along cable
	mvecs.clear();
	double y;
	double x;
	double z;
	nb_nodes = nb_elems * 2 + 1;
	double ds = length / ((double)nb_nodes - 1);
	for (int i = 0; i < nb_nodes; ++i) {
		x = x0 + a*asinh((L0 + ds*i) / a);
		y = 0.001 + a*(cosh(x / a) - 1);
		z = 0;
		//x = L0 + ds*i;
		//y = 1.5+0.00001*x;
		mvecs.push_back(ChVector<>(x, y, z));
		elems_length.push_back(ds * 2);
	}
}

void cppCable::buildNodes() {
	nodes.clear();
  nb_nodes = 2*nb_elems+1;
	ChVector<> dir;  // direction of node
	ChVector<> norm1; // normal of node direction
	ChVector<> norm2; // normal of node direction and norm1
	ChCoordsys<> coordsys;  // coordinate system
	ChVector<> ref = ChVector<>(1., 0., 0.);
	std::shared_ptr<ChNodeFEAxyzDD> node;
	// first node
	dir = mvecs[1] - mvecs[0];
	dir.Normalize();
	//plane = dir.x()*(x - mvecs[0].x()) + dir.y()*(y - mvecs[0].y()) + dir.z()*(z - mvecs[0].z());
	if (dir.x() == 1 && dir.y() == 0 && dir.z() == 0) {
		ref = ChVector<>(0., 0., -1.);
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
	for (int i = 1; i < nb_nodes - 1; ++i) {
		dir = mvecs[i + 1] - mvecs[i - 1];
		dir.Normalize();
		if (dir.x() == 1 && dir.y() == 0 && dir.z() == 0) {
			ref = ChVector<>(0., 0., -1.);
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
	if (dir.x() == 1 && dir.y() == 0 && dir.z() == 0) {
		ref = ChVector<>(0., 0., -1.);
		ref.Normalize();
	}
	else {
		ref = ChVector<>(1., 0., 0.);
	}
	norm1 = dir % ref;
	norm1.Normalize();
	norm2 = dir % norm1;
	norm2.Normalize();
	node = std::make_shared<ChNodeFEAxyzDD>(mvecs[nb_nodes - 1], norm1, norm2);
	nodes.push_back(node);
}

void cppCable::buildElements() {
	// auto loadcontainer = std::make_shared<ChLoadContainer>();
	// system.Add(loadcontainer);
	// build elements
	elems.clear();
  nb_nodes = nb_elems*2+1;
	for (int i = 1; i < nb_nodes - 1; ++i) {
		if (i % 2 != 0) {
			auto element = std::make_shared<ChElementBeamANCF>();
			elems.push_back(element);
			element->SetNodes(nodes[i - 1], nodes[i + 1], nodes[i]);
      double elem_length = (nodes[i]->GetPos()-nodes[i-1]->GetPos()).Length()+(nodes[i+1]->GetPos()-nodes[i]->GetPos()).Length();
			element->SetDimensions(elem_length, d, d);
			element->SetMaterial(mmaterial_cable);
			element->SetAlphaDamp(0.0004);
			element->SetGravityOn(true);
			element->SetStrainFormulation(ChElementBeamANCF::StrainFormulation::CMNoPoisson);
		}
	}
}

void cppCable::buildMesh() {
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

void cppCable::setFluidVelocityAtNodes(std::vector<ChVector<>> vel) {
	fluid_velocity = vel;
}

void cppCable::setDragForce() {
	ChVector<> u_ch;  // velocity from chrono
	ChVector<> u_prot;  // velocity from proteus
	ChVector<> u_rel;  // relative velocity of node with surrounding fluid
	ChVector<> t_dir;  // tangent at node
	ChVector<> Fd_trans;  // transversal drag force
	ChVector<> Fd_tang;  // trangential drag force
	ChVector<> Fd;  // total drag force
	std::shared_ptr<ChVector<>> force_drag;
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
		u_rel = u_ch - u_prot;
		// CAREFUL HERE: ChBeamElementANCF, GetD() does not give direction but normal
		t_dir = nodes[i]->GetD();
		// transverse drag coefficient
		double C_dn = 0.;  // transversal drag coeff
		double C_dt = 0.;  // tangential drag coeff
		Fd_trans = 0.5*rho*C_dn*length_elem*((u_rel^t_dir)*t_dir - u_rel).Length()*((u_rel^t_dir)*t_dir - u_rel);
		Fd_tang = 0.5*rho*C_dt*M_PI*length_elem*((-u_rel^t_dir)*t_dir).Length()*((-u_rel^t_dir)*t_dir);
		Fd = Fd_trans + Fd_tang;
		forces_drag.push_back(std::make_shared<ChVector<>>(Fd.x(), Fd.y(), Fd.z()));
	}
}

void cppCable::getBuoyancyForce() {
	/*
	Here probably use the SetDensity() attributes of element to set the new density
	TO CHANGE !!!
	*/
	//double vof;  // get vof (from proteus)
	//double rho_w;  // density of water (from proteus)
	//double rho_a;  // density of air (from proteus)
	//double rho_elem;  // density of elem
	//for (int i = 0; i < elems.size(); ++i) {
	//	rho_elem = rho - (rho_a*vof + rho_w*(1. - vof));
	//	elems[i]->GetSection()->SetDensity(rho_elem);
	//}
}

cppMultiSegmentedCable * newMoorings(
	ChSystemDEM& system,
	std::shared_ptr<ChMesh> mesh,
	std::vector<double> length,
	std::vector<int> nb_elems,
	std::vector<double> d,
	std::vector<double> rho,
	std::vector<double> E)
{
	return new cppMultiSegmentedCable(system, mesh, length, nb_elems, d, rho, E);
}

cppMesh * newMesh(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh) {
	return new cppMesh(system, mesh);
}







/* class cppSurfaceBoxNodesCloud { */
/*   ChSystemDEM& system; */
/*   std::shared_ptr<ChMesh>; */
/*   ChVector<> position; */
/*   ChVector<> dimensions; */
/*   std::shared_ptr<ChBodyEasyBox> box; */
/*   std::shared_ptr<ChMaterialSurfaceDEM> material; */
/*   std::shared_ptr<ChContactSurfaceNodeCloud> contact_cloud; */
/*   cppSurfaceBoxNodesCloud(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, Chvector position, ChVector dimensions); */
/*   void setNodesSize(double size); */
/* } */

/* cppSurfaceBoxNodesCloud::cppSurfaceBoxNodesCloud(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, Chvector position, ChVector dimensions) { */
/* 	// Create a surface material to be shared with some objects */
/*   system = system; */
/* 	material = std::make_shared<ChMaterialSurfaceDEM>(); */
/* 	material->SetYoungModulus(2e4); */
/* 	material->SetFriction(0.3f); */
/* 	material->SetRestitution(0.2f); */
/* 	material->SetAdhesion(0); */
/* 	contact_cloud = std::make_shared<ChContactSurfaceNodeCloud>(); */
/*   mesh = mesh; */
/* 	mesh->AddContactSurface(contact_cloud); */
/* 	// Must use this to 'populate' the contact surface.  */
/* 	// Use larger point size to match beam section radius */
/* 	contact_cloud->AddAllNodes(0.001); */

/* 	// Use our DEM surface material properties  */
/* 	contact_cloud->SetMaterialSurface(material); */

/* 	box = std::make_shared<ChBodyEasyBox>( */
/*     dimensions.x(), dimensions.y(), dimensions.z(),  // x,y,z size */
/* 		1000,       // density */
/* 		true,       // visible */
/* 		true        // collide */
/* 		); */

/* 	system.Add(box); */

/* 	box->SetBodyFixed(true); */
/* 	box->SetPos(position); */

/* 	// Use our DEM surface material properties  */
/* 	box->SetMaterialSurface(mysurfmaterial); */
/* } */

/* void SurfaceBoxNodesCloud::setNodesSize(double size) { */
/*   contact_cloud->AddAllNodes(size); */
/* } */

/* cppSurfaceBoxNodesCloud * newSurfaceBoxNodesCloud(ChSystemDEM& system, std::shared_ptr<ChMesh> mesh, ChVector<> position, ChVector<> dimensions) { */
/* 	return new cppSurfaceBoxNodesCloud(system, mesh, position, dimensions); */
/* } */
