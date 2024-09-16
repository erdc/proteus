//#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <string>
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactMaterial.h"
#include "chrono/physics/ChContactMaterialSMC.h"
#include "chrono/fea/ChElementBeamANCF_3333.h"
#include "chrono/fea/ChElementCableANCF.h"
#include "chrono/fea/ChElementBeamEuler.h"
#include "chrono/fea/ChBeamSection.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChLinkNodeNode.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChLinkNodeSlopeFrame.h"
#include "chrono/fea/ChLoadsBeam.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverPMINRES.h"


//using namespace std;
using namespace chrono;
using namespace fea;


// override some functions of ChElement

class ChElementCableANCFmod : public ChElementCableANCF {
private:
  ChVectorN<double, 12> m_GenForceVec0;
  virtual void SetupInitial(ChSystem* system) override {
    assert(GetSection());

    // Compute rest length, mass:
    //double length2 = (nodes[1]->GetX0() - nodes[0]->GetX0()).Length();
    //this->mass = this->length * GetSection()->Area * GetSection()->density;
    this->mass = this->length * GetDensity();

    // Here we calculate the internal forces in the initial configuration
    // Contribution of initial configuration in elastic forces is automatically subtracted
    ChVectorDynamic<> FVector0(12);
    FVector0.setZero();
    m_GenForceVec0.setZero();
    ComputeInternalForces(FVector0);
    m_GenForceVec0 = FVector0;

    // Compute mass matrix
    ComputeMassMatrix();
  };
};

class ChElementBeamEulermod : public ChElementBeamEuler {
private:
  ChQuaternion<> q_element_ref_rot;
  virtual void SetupInitial(ChSystem* system) override {
    assert(GetSection());

    // Compute rest length, mass:
    //this->length = (nodes[1]->GetX0().GetPos() - nodes[0]->GetX0().GetPos()).Length();
    //this->mass = this->length * GetSection()->Area * GetSection()->density;
    this->mass = this->length * GetDensity();

    // Compute initial rotation
    ChMatrix33<> A0;
    auto node0 = GetNodeA();
    auto node1 = GetNodeB();
    ChVector3d mXele = node1->GetX0().GetPos() - node0->GetX0().GetPos();
    ChVector3d myele = node0->GetX0().GetRotMat().GetAxisY();
    A0.SetFromAxisX(mXele, myele);
    q_element_ref_rot = A0.GetQuaternion();

    // Compute local stiffness matrix:
    ComputeStiffnessMatrix();
  };
};


class MyLoaderTriangular : public ChLoaderUdistributed {
public:
  // Useful: a constructor that also sets ChLoadable
  ChVector3d Fa;
  ChVector3d Fb;
  MyLoaderTriangular(std::shared_ptr<ChLoadableU> mloadable):
    ChLoaderUdistributed(mloadable) {
    Fa = ChVector3d(0.,0.,0.);
    Fb = ChVector3d(0.,0.,0.);
  };
  // Compute F=F(u)
  // This is the function that you have to implement. It should return the
  // load at U. For Eulero beams, loads are expected as 6-rows vectors, containing
  // a wrench: forceX, forceY, forceZ, torqueX, torqueY, torqueZ.
  virtual void ComputeF(const double U,
                        ChVectorDynamic<>& F,
                        ChVectorDynamic<>* state_x,
                        ChVectorDynamic<>* state_w
                        ) {
    double Fy_max = 0.005;
    ChVector3d force = Fa*abs(-1+U)/2.+Fb*(1+U)/2.;
    F(0) = force.x();
    F(1) = force.y();
    F(2) = force.z();
  }
  void SetF(ChVector3d Fa_in, ChVector3d Fb_in) {
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
  std::shared_ptr<ChSystem> system;  // global system
  std::shared_ptr<ChMesh> mesh;  // mesh
  int nb_elems;   // number of nodes along cable
  std::vector<double> length_per_elem;  // length of each element on the cable
  double Cd_axial;  // drag coeff in axial direction
  double Cd_normal;  // drag coeff in normal direction
  double Cm_axial;  // added mass coeff in axial direction
  double Cm_normal;  // added mass coeff in normal direction
  std::vector<ChVector3d> mvecs;  // vectors (nodes coordinates)
  std::vector<ChVector3d> mvecs_tangents;  // vectors (tangents at nodes coordinates)
  std::vector<ChVector3d> mvecs_middle;
  std::vector<ChVector3d> mdirs;  // vectors (nodes coordinates)
  double d, rho, E, length;  // diameter, density, Young's modulus, length of cable
  double A0; // unstretched diameter
  double L0 = 0;  // initial length along cable
  double Iyy;
  int nb_nodes;
  bool applyDrag;
  bool applyAddedMass;
  bool applyBuoyancy;
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
  std::vector<ChVector3d> fluid_velocity;
  std::vector<ChVector3d> fluid_acceleration;
  std::vector<double> fluid_density;
  std::vector<double> nodes_density; // density of (cable-fluid) at nodes
  cppCable(std::shared_ptr<ChSystem> system, std::shared_ptr<ChMesh> mesh, double length,
           int nb_elems, double d, double rho, double E, double L0, std::string beam_type);  // constructor
  void setFluidVelocityAtNodes(std::vector<ChVector3d> vel);
  void setFluidAccelerationAtNodes(std::vector<ChVector3d> acc);
  void setFluidDensityAtNodes(std::vector<double> vof);
  std::vector<std::shared_ptr<ChVector3d>> getNodalPositions();
  std::vector<std::shared_ptr<ChVector3d>> forces_drag;
  std::vector<std::shared_ptr<ChVector3d>> forces_addedmass;
  /* std::vector<std::shared_ptr<ChLoadBeamWrenchDistributed>> elems_loads_distributed; */
  std::vector<std::shared_ptr<ChLoad>> elems_loads_triangular;
  std::vector<std::shared_ptr<ChLoad>> elems_loads_volumetric;
  /* std::vector<std::shared_ptr<ChLoadBeamWrench>> elems_loads; */
  void buildVectors();  // builds location vectors for the nodes
  void buildNodes(bool last_node);  // builds the nodes for the mesh
  void buildNodesBeamEuler(bool last_node);  // builds the nodes for the mesh
  void buildNodesCableANCF(bool last_node);  // builds the nodes for the mesh
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
  void setRestLengthPerElement(std::vector<double> length_array);
  void setIyy(double Iyy_in);
};

class cppMultiSegmentedCable {
public:
  std::shared_ptr<ChSystem> system;  // global system
  std::string beam_type;
  std::shared_ptr<ChContactMaterialSMC> mysurfmaterial;
  std::shared_ptr<ChMesh> mesh;  // mesh
  std::shared_ptr<ChBody> fairleadd;
  std::shared_ptr<ChLinkNodeFrame> fairlead2;
  std::vector<int> nb_nodes;   // number of nodes along cable
  std::vector<int> nb_elems;   // number of nodes along cable
  std::vector<ChVector3d> mvecs;  // vectors (nodes coordinates)
  std::vector<std::shared_ptr<cppCable>> cables;
  std::shared_ptr<ChContactMaterialSMC> contact_material;  // mesh
  std::vector<double>  d;
  std::vector<double> rho;
  std::vector<double> E;
  std::vector<double> length;  // diameter, density, Young's modulus, length of cable
  std::vector<std::shared_ptr<ChNodeFEAxyzD>> nodes;  // array nodes coordinates and direction
  std::vector<std::shared_ptr<ChNodeFEAxyzDD>> nodesDD;  // array nodes coordinates and direction
  std::vector<std::shared_ptr<ChNodeFEAxyzrot>> nodesRot;  // array nodes coordinates and direction
  std::vector<ChVector3d> fluid_velocity;
  std::vector<ChVector3d> fluid_acceleration;
  std::vector<double> fluid_density;
  std::vector<std::shared_ptr<ChElementCableANCF>> elemsCableANCF;  // array of elements */
  std::vector<std::shared_ptr<ChElementBeamEuler>> elemsBeamEuler;  // array of elements */
  std::shared_ptr<ChLinkBase> constraint_front;
  std::shared_ptr<ChLinkBase> constraint_back;
  std::vector<std::shared_ptr<ChVector3d>> forces_drag;
  std::vector<std::shared_ptr<ChVector3d>> forces_addedmass;
  std::shared_ptr<ChBody> body_back;
  std::shared_ptr<ChBody> body_front;
  int nb_nodes_tot;
  int nb_elems_tot;
  bool nodes_built;
  bool elems_built;
  bool nodes_chlink;
  cppMultiSegmentedCable(std::shared_ptr<ChSystem> system,
                         std::shared_ptr<ChMesh> mesh,
                         std::vector<double> length,
                         std::vector<int> nb_nodes,
                         std::vector<double> d,
                         std::vector<double> rho,
                         std::vector<double> E,
                         std::string beam_type);
  void setFluidVelocityAtNodes(std::vector<ChVector3d> vel);
  void setFluidAccelerationAtNodes(std::vector<ChVector3d> vel);
  void setFluidDensityAtNodes(std::vector<double> dens);
  void updateDragForces();
  void updateAddedMassForces();
  void applyForces();
  std::vector<std::shared_ptr<ChVector3d>> getNodalPositions();
  void buildNodes();
  void buildElements();
  void buildCable();  // builds the multi-segmented cable
  void getForceFairlead();

  void attachBackNodeToBody(std::shared_ptr<ChBody> body);
  void attachFrontNodeToBody(std::shared_ptr<ChBody> body);
  void setContactMaterial(std::shared_ptr<ChContactMaterialSMC> material);
  void buildNodesCloud();
  ChVector3d getTensionElement(int i, double eta);
};


cppMultiSegmentedCable::cppMultiSegmentedCable(std::shared_ptr<ChSystem> system,
                                               std::shared_ptr<ChMesh> mesh,
                                               std::vector<double> length,
                                               std::vector<int> nb_elems,
                                               std::vector<double> d,
                                               std::vector<double> rho,
                                               std::vector<double> E,
                                               std::string beam_type="CableANCF"):
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
  elems_built = false;
  nodes_chlink = true;  // build links (true) or link elements directly (false)

  std::shared_ptr<cppCable> segment;
  double L0 = 0;
  for (int i = 0; i < length.size(); ++i) {
    segment = std::make_shared<cppCable>(system,
                                         mesh,
                                         length[i],
                                         nb_elems[i],
                                         d[i],
                                         rho[i],
                                         E[i],
                                         L0,
                                         beam_type);
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
      nodesRot.insert(nodesRot.end(),
                      cables[i]->nodesRot.begin(),
                      cables[i]->nodesRot.end());
      nb_nodes_tot = nodesRot.size();
    }
    else if (beam_type == "CableANCF") {
      if (nodes_chlink == false && i < cables.size()-1) {
        cables[i]->buildNodes(false);
      }
      else {
        cables[i]->buildNodes(true);
      }
      nodes.insert(nodes.end(),
                   cables[i]->nodes.begin(),
                   cables[i]->nodes.end());
      nb_nodes_tot = nodes.size();
    }
  }
  nodes_built = true;
}

void cppMultiSegmentedCable::buildElements() {
  elemsCableANCF.clear();
  elemsBeamEuler.clear();
  for (int i = 0; i < cables.size(); ++i) {
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
      nb_elems_tot = elemsCableANCF.size();
    }
    else if (beam_type == "BeamEuler") {
      elemsBeamEuler.insert(elemsBeamEuler.end(), cables[i]->elemsBeamEuler.begin(), cables[i]->elemsBeamEuler.end());
      nb_elems_tot = elemsBeamEuler.size();
    }
  }
  elems_built = true;
}

void cppMultiSegmentedCable::buildCable() {
  /* builds all cable segments and updates their link
     (no duplicate node added to mesh) */
  if (nodes_built == false) {
    buildNodes();
  }
  if (elems_built == false) {
    buildElements();
  }
  for (int i = 0; i < cables.size(); ++i) {
    if (i < cables.size()-1 && nodes_chlink == false) {
      cables[i]->buildMesh(false);
    }
    else {
      cables[i]->buildMesh(true);
    }
    if (nodes_chlink == true) {
      if (i>0) {
        if (beam_type == "BeamEuler") {
          auto con1 = chrono_types::make_shared<ChLinkMateSpherical>();
          auto nodeA = cables[i]->nodesRot.front();
          auto nodeB = cables[i-1]->nodesRot.back();
          con1->Initialize(nodeA, nodeB, false, nodeA->GetPos(), nodeA->GetPos());
          system->Add(con1);
        }
        else if (beam_type == "CableANCF") {
          auto con1 = chrono_types::make_shared<ChLinkNodeNode>();
          auto nodeA = cables[i]->nodes.front();
          auto nodeB = cables[i-1]->nodes.back();
          con1->Initialize(nodeA, nodeB);
          system->Add(con1);
        }
      }
    }
    forces_drag.insert(forces_drag.end(), cables[i]->forces_drag.begin(), cables[i]->forces_drag.end());
    forces_addedmass.insert(forces_addedmass.end(), cables[i]->forces_addedmass.begin(), cables[i]->forces_addedmass.end());
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
    fluid_velocity.push_back(ChVector3d(0.,0.,0.));
    fluid_acceleration.push_back(ChVector3d(0.,0.,0.));
    fluid_density.push_back(0.);
  }
  setFluidVelocityAtNodes(fluid_velocity);
  setFluidAccelerationAtNodes(fluid_acceleration);
  setFluidDensityAtNodes(fluid_density);
}

void cppMultiSegmentedCable::setFluidAccelerationAtNodes(std::vector<ChVector3d> acc) {
  fluid_acceleration = acc;
  int node_nb = 0;
  int node_nb_prev = node_nb;
  for (int i = 0; i < cables.size(); ++i) {
    if (beam_type == "BeamEuler") {
      node_nb += cables[i]->nodesRot.size();
    }
    else if (beam_type == "CableANCF") {
      node_nb += cables[i]->nodes.size();
    }
    std::vector<ChVector3d> fluid_acc(fluid_acceleration.begin()+node_nb_prev,
                                      fluid_acceleration.begin()+node_nb);
    cables[i]->setFluidAccelerationAtNodes(fluid_acc);
    node_nb_prev = node_nb;
  }
}

void cppMultiSegmentedCable::setFluidVelocityAtNodes(std::vector<ChVector3d> vel) {
  fluid_velocity = vel;
  int node_nb = 0;
  int node_nb_prev = node_nb;
  for (int i = 0; i < cables.size(); ++i) {
    if (beam_type == "BeamEuler") {
      node_nb += cables[i]->nodesRot.size();
    }
    else if (beam_type == "CableANCF") {
      node_nb += cables[i]->nodes.size();
    }
    std::vector<ChVector3d> fluid_vel(fluid_velocity.begin()+node_nb_prev,
                                      fluid_velocity.begin()+node_nb);
    cables[i]->setFluidVelocityAtNodes(fluid_vel);
    node_nb_prev = node_nb;
  }
}

void cppMultiSegmentedCable::setFluidDensityAtNodes(std::vector<double> dens) {
  fluid_density = dens;
  int node_nb = 0;
  int node_nb_prev = 0;
  for (int i = 0; i < cables.size(); ++i) {
    if (beam_type == "BeamEuler") {
      node_nb += cables[i]->nodesRot.size();
    }
    else if (beam_type == "CableANCF") {
      node_nb += cables[i]->nodes.size();
    }
    std::vector<double> fluid_dens(fluid_density.begin()+node_nb_prev,
                                   fluid_density.begin() + node_nb);
    cables[i]->setFluidDensityAtNodes(fluid_dens);
    node_nb_prev = node_nb;
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

ChVector3d cppMultiSegmentedCable::getTensionElement(int i, const double eta=0.) {

  auto force = ChVector3d();
  auto torque = ChVector3d();
  if (beam_type == "CableANCF") {
    elemsCableANCF[i]->EvaluateSectionForceTorque(eta,
                                                  force,
                                                  torque);
    elemsCableANCF[i]->EvaluateSectionStrain(eta,
                                             force);
  }
  else if (beam_type == "BeamEuler") {
    auto mat2 = ChMatrixDynamic<>();
    elemsBeamEuler[i]->EvaluateSectionForceTorque(eta,
                                                  force,
                                                  torque);
  }
  return force;
}

std::vector<std::shared_ptr<ChVector3d>> cppMultiSegmentedCable::getNodalPositions() {
  std::vector<std::shared_ptr<ChVector3d>> nodal_positions;
  for (int i = 0; i < nodes.size(); ++i) {
    auto pos = nodes[i]->GetPos();
    double x = pos.x();
    double y = pos.y();
    double z = pos.z();
    auto nodal_position = chrono_types::make_shared<ChVector3d>(x, y, z);
    nodal_positions.push_back(nodal_position);
  }
  return nodal_positions;
}

void cppMultiSegmentedCable::attachBackNodeToBody(std::shared_ptr<ChBody> body) {
  if (beam_type == "BeamEuler") {
    auto constraint = chrono_types::make_shared<ChLinkMateSpherical>();
    constraint->Initialize(nodesRot.back(), body, false, nodesRot.back()->GetPos(), nodesRot.back()->GetPos());
    system->Add(constraint);
    body_back = body;
    constraint_back = constraint;
  }
  else {
    auto constraint = chrono_types::make_shared<ChLinkNodeFrame>();
    constraint->Initialize(nodes.back(), body);
    system->Add(constraint);
    body_back = body;
    constraint_back = constraint;
  }
};

void cppMultiSegmentedCable::attachFrontNodeToBody(std::shared_ptr<ChBody> body) {
  if (beam_type == "BeamEuler") {
    auto constraint = chrono_types::make_shared<ChLinkMateSpherical>();
    constraint->Initialize(nodesRot.front(), body, false, nodesRot.front()->GetPos(), nodesRot.front()->GetPos());
    system->Add(constraint);
    body_front = body;
    constraint_front = constraint;
  }
  else if (beam_type == "CableANCF") {
    auto constraint = chrono_types::make_shared<ChLinkNodeFrame>();
    constraint->Initialize(nodes.front(), body);
    system->Add(constraint);
    body_front = body;
    constraint_front = constraint;
  }
};

void cppMultiSegmentedCable::setContactMaterial(std::shared_ptr<ChContactMaterialSMC> material) {
  contact_material = material;
};

void cppMultiSegmentedCable::buildNodesCloud() {
  if (contact_material) {
    // Use DEM surface material properties
    auto contact_cloud = chrono_types::make_shared<ChContactSurfaceNodeCloud>(contact_material);
    //auto contact_cloud = chrono_types::make_shared<ChContactSurfaceNodeCloud>();
    mesh->AddContactSurface(contact_cloud);
    //contact_cloud->SetMaterialSurface(contact_material);
    // add cable nodes to cloud
    for (int i = 0; i < cables.size(); ++i) {
      cables[i]->addNodestoContactCloud(contact_cloud);
    }
  }
};


cppCable::cppCable(std::shared_ptr<ChSystem> system, // system in which the cable belong
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
  /* Iyy = 1e-12; */
  applyDrag = true;
  applyAddedMass = true;
  applyBuoyancy = true;
  for (int i = 0; i < nb_elems; i++) {
    length_per_elem.push_back(length/nb_elems);
  }
  if (beam_type == "CableANCF") {
    msection_cable = chrono_types::make_shared<ChBeamSectionCable>();
    msection_cable->SetDiameter(d);
    msection_cable->SetYoungModulus(E);
    msection_cable->SetDensity(rho);
    /* msection_cable->SetInertia(Iyy); */
    Iyy = msection_cable->GetInertia();
  }
  else if (beam_type == "BeamEuler") {
    msection_advanced = chrono_types::make_shared<ChBeamSectionAdvanced>();
    msection_advanced->SetYoungModulus(E);
    msection_advanced->SetShearModulus(1e-6);
    msection_advanced->SetDensity(rho);
    msection_advanced->SetAsCircularSection(d);
    /* msection_advanced->SetIyy(Iyy); */
    /* msection_advanced->SetIzz(Iyy); */
    Iyy = msection_advanced->GetIyy();
  }
}

void cppCable::setIyy(double Iyy_in) {
  Iyy = Iyy_in;
  if (beam_type == "CableANCF") {
    msection_cable->SetInertia(Iyy);
  }
  else if (beam_type == "BeamEuler") {
    msection_advanced->SetIyy(Iyy);
    msection_advanced->SetIzz(Iyy);
  }
}

void cppCable::buildNodes(bool last_node=true) {
  if (beam_type == "CableANCF") {buildNodesCableANCF(last_node);}
  else if (beam_type == "BeamEuler") {buildNodesBeamEuler(last_node);}
}

void cppCable::buildNodesBeamEuler(bool last_node) {
  nodesRot.clear();
  forces_drag.clear();
  forces_addedmass.clear();
  std::shared_ptr<ChNodeFEAxyzrot> node;
  ChVector3d dir;  // direction of node
  ChVector3d ref = ChVector3d(1.,0.,0.);
  ChQuaternion<> frame_quat;
  for (int i = 0; i < mvecs.size() - 1; ++i) {
    dir = mvecs_tangents[i];
    dir.Normalize();
    double ang = acos(dir^ref);  // inner product
    auto axis = ref%dir; // cross product
    frame_quat.SetFromAngleAxis(ang, axis);
    node = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(mvecs[i],
								frame_quat));
    nodesRot.push_back(node);
    std::shared_ptr<ChVector3d> drag0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    std::shared_ptr<ChVector3d> am0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    forces_drag.push_back(drag0);
    forces_addedmass.push_back(am0);
  }  // last node
  if (last_node == true) {
    dir = mvecs_tangents[mvecs.size()-1];
    dir.Normalize();
    double ang = -acos(dir^ref);  // inner product
    auto axis = ref%dir; // cross product
    frame_quat.SetFromAngleAxis(ang, axis);
    node = chrono_types::make_shared<ChNodeFEAxyzrot>(ChFrame<>(mvecs[mvecs.size()-1],
								frame_quat));
    nodesRot.push_back(node);
    nb_nodes = nodesRot.size();
    nb_elems = nb_nodes-1;
    std::shared_ptr<ChVector3d> drag0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    std::shared_ptr<ChVector3d> am0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    forces_drag.push_back(drag0);
    forces_addedmass.push_back(am0);
  }
  else {
    nb_nodes = nodesRot.size();
    nb_elems = nb_nodes;
  }
}

void cppCable::buildNodesCableANCF(bool last_node) {
  nodes.clear();
  forces_drag.clear();
  forces_addedmass.clear();
  std::shared_ptr<ChNodeFEAxyzD> node;
  ChVector3d dir;  // direction of node
  ChCoordsys<> coordsys;  // coordinate system of node
  for (int i = 0; i < mvecs.size() - 1; ++i) {
    dir = mvecs_tangents[i];
    dir.Normalize();
    node = chrono_types::make_shared<ChNodeFEAxyzD>(mvecs[i], dir);
    nodes.push_back(node);
    std::shared_ptr<ChVector3d> drag0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    std::shared_ptr<ChVector3d> am0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    forces_drag.push_back(drag0);
    forces_addedmass.push_back(am0);
  }  // last node
  if (last_node == true) {
    dir = mvecs_tangents[mvecs_tangents.size()-1];
    dir.Normalize();
    node = chrono_types::make_shared<ChNodeFEAxyzD>(mvecs[mvecs.size()-1], dir);
    nodes.push_back(node);
    nb_nodes = nodes.size();
    nb_elems = nb_nodes-1;
    std::shared_ptr<ChVector3d> drag0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    std::shared_ptr<ChVector3d> am0 = chrono_types::make_shared<ChVector3d>(0.,0.,0.);
    forces_drag.push_back(drag0);
    forces_addedmass.push_back(am0);
  }
  else {
    nb_nodes = nodes.size();
    nb_elems = nb_nodes;
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
  auto loadcontainer = chrono_types::make_shared<ChLoadContainer>();
  system->Add(loadcontainer);
  // build elements
  elemsCableANCF.clear();
  /* elems_loads_distributed.clear(); */
  elems_loads_volumetric.clear();
  elems_loads_triangular.clear();
  /* elems_loads.clear(); */
  for (int i = 0; i < nb_elems; ++i) {
    auto element = chrono_types::make_shared<ChElementCableANCFmod>();
    auto load_distributed = chrono_types::make_shared<ChLoadBeamWrenchDistributed>(element);
    auto load = chrono_types::make_shared<ChLoadBeamWrench>(element);
    auto gravity = chrono_types::make_shared<ChLoaderGravity>(element);
    gravity->SetGravitationalAcceleration(system->GetGravitationalAcceleration());
    std::shared_ptr<ChLoad> loadtri(new ChLoad(gravity));
    auto load_volumetric = chrono_types::make_shared<ChLoad>(gravity);
    /* loadcontainer->Add(load_distributed); */
    /* loadcontainer->Add(load); */
    loadcontainer->Add(loadtri);  // do not forget to add the load to the load container.
    loadcontainer->Add(load_volumetric);
    elemsCableANCF.push_back(element);
    /* elems_loads_distributed.push_back(load_distributed); */
    /* elems_loads.push_back(load); */
    elems_loads_triangular.push_back(loadtri);
    elems_loads_volumetric.push_back(load_volumetric);
    element->SetSection(msection_cable);
    element->SetRestLength(length_per_elem[i]);
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
  auto loadcontainer = chrono_types::make_shared<ChLoadContainer>();
  system->Add(loadcontainer);
  // build elements
  elemsBeamEuler.clear();
  /* elems_loads_distributed.clear(); */
  elems_loads_triangular.clear();
  /* elems_loads.clear(); */
  for (int i = 0; i < nodesRot.size() - 1; ++i) {
    auto element = chrono_types::make_shared<ChElementBeamEulermod>();
    auto load_distributed = chrono_types::make_shared<ChLoadBeamWrenchDistributed>(element);
    auto load = chrono_types::make_shared<ChLoadBeamWrench>(element);
    auto gravity = chrono_types::make_shared<ChLoaderGravity>(element);
    gravity->SetGravitationalAcceleration(system->GetGravitationalAcceleration());
    std::shared_ptr<ChLoad> loadtri(new ChLoad(gravity));
    auto load_volumetric = chrono_types::make_shared<ChLoad>(gravity);
    /* loadcontainer->Add(load_distributed); */
    /* loadcontainer->Add(load); */
    loadcontainer->Add(loadtri);  // do not forget to add the load to the load container.
    loadcontainer->Add(load_volumetric);
    elemsBeamEuler.push_back(element);
    /* elems_loads_distributed.push_back(load_distributed); */
    /* elems_loads.push_back(load); */
    elems_loads_triangular.push_back(loadtri);
    elems_loads_volumetric.push_back(load_volumetric);
    element->SetSection(msection_advanced);
    element->SetRestLength(length/nb_elems);
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
void cppCable::setFluidAccelerationAtNodes(std::vector<ChVector3d> acc) {
  fluid_acceleration = acc;
}

void cppCable::setFluidVelocityAtNodes(std::vector<ChVector3d> vel) {
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

void cppCable::setRestLengthPerElement(std::vector<double> length_array) {
  length_per_elem = length_array;
}

void cppCable::setDragForce() {
  /*
   * setFluidVelocityAtNodes and setFluidDensityAtNodes
   * must be called before this function
   */
  ChVector3d u_ch;  // velocity from chrono
  ChVector3d u_prot;  // velocity from proteus
  ChVector3d u_rel;  // relative velocity of node with surrounding fluid
  ChVector3d t_dir;  // tangent at node
  ChVector3d Fd_a;  // axial (tangential) drag force
  ChVector3d Fd_n;  // normal(transversal) drag force
  ChVector3d Fd;  // total drag force
  ChVector3d Va;
  ChVector3d Vn;
  double rho_f;
  // clear current drag forces
  double length_elem = length / (nb_nodes - 1);
  for (int i = 0; i < nb_nodes; ++i) {
    if (beam_type == "CableANCF") {
      t_dir = nodes[i]->GetSlope1();
      u_ch = nodes[i]->GetPosDt();
    }
    else if (beam_type == "BeamEuler") {
      t_dir = nodesRot[i]->GetRot().GetVector();
      u_ch = nodesRot[i]->GetPosDt();
    }
    // get velocity u_prot from proteus // TO CHANGE !!
    double ux_prot = fluid_velocity[i][0];
    double uy_prot = fluid_velocity[i][1];
    double uz_prot = fluid_velocity[i][2];
    u_prot = ChVector3d(ux_prot, uy_prot, uz_prot);
    u_rel = u_prot - u_ch;
    // CAREFUL HERE: ChBeamElementANCF, GetD() does not give direction but normal
    rho_f = fluid_density[i];
    double dot = u_rel^t_dir;
    Va = t_dir*dot;
    Vn = u_rel-Va;
    Fd_a = 0.5*rho_f*Cd_axial*M_PI*d*Va.Length()*Va;//(force per unit length)
    Fd_n = 0.5*rho_f*Cd_normal*d*Vn.Length()*Vn;//(force per unit length)
    Fd = Fd_a + Fd_n;
    forces_drag[i]->Set(Fd);
  }
}


void cppCable::setAddedMassForce() {
  /*
   * setFluidVelocityAtNodes and setFluidDensityAtNodes
   * must be called before this function
   */
  ChVector3d a_ch;  // acceleration from chrono
  ChVector3d a_prot;  // acceleration from proteus
  ChVector3d a_rel;  // relative acceleration of node with surrounding fluid
  ChVector3d t_dir;  // tangent at node
  ChVector3d Fm_a;  // axial (tangential) added mass force
  ChVector3d Fm_n;  // normal(transversal) added mass force
  ChVector3d Fm_f;  // fluid part added mass force
  ChVector3d Fm;  // total added mass force
  ChVector3d Va;
  ChVector3d Vn;
  double rho_f;
  // clear current drag forces
  double length_elem = length / (nb_nodes - 1);
  for (int i = 0; i < nb_nodes; ++i) {
    if (beam_type == "CableANCF") {
      t_dir = nodes[i]->GetSlope1();
      a_ch = nodes[i]->GetPosDt2();
    }
    else if (beam_type == "BeamEuler") {
      t_dir = nodesRot[i]->GetRot().GetVector();
      a_ch = nodesRot[i]->GetPosDt2();
    }
    // get velocity u_prot from proteus // TO CHANGE !!
    double ax_prot = fluid_acceleration[i][0];
    double ay_prot = fluid_acceleration[i][1];
    double az_prot = fluid_acceleration[i][2];
    a_prot = ChVector3d(ax_prot, ay_prot, az_prot);
    a_rel = a_prot - a_ch;
    rho_f = fluid_density[i];
    double dot = a_rel^t_dir;
    Va = t_dir*dot;
    Vn = a_rel-Va;
    Fm_a = rho_f*Cm_axial*M_PI*d*d/4.*Va;//(force per unit length)
    Fm_n = rho_f*Cm_normal*M_PI*d*d/4.*Vn;//(force per unit length)
    Fm_f = rho_f*M_PI*d*d/4.*a_prot;
    Fm = Fm_a + Fm_n + Fm_f;
    forces_addedmass[i]->Set(Fm);
  }
}

void cppCable::applyForces() {
  for (int i = 0; i < nb_nodes-1; ++i) {
    ChVector3d Fa = ChVector3d(0.,0.,0.);
    ChVector3d Fb = ChVector3d(0.,0.,0.);
    if (applyDrag == true) {
      Fa = Fa+*forces_drag[i].get();
      Fb = Fb+*forces_drag[i+1].get();
    }
    if (applyAddedMass == true) {
      Fa = Fa+*forces_addedmass[i].get();
      Fb = Fb+*forces_addedmass[i+1].get();
    }
    /* elems_loads_triangular[i]->loader.SetF(Fa, Fb); */
    // buoyancy
    if (applyBuoyancy == true) {
      if (mesh->GetAutomaticGravity() == true) {
        /* elems_loads_volumetric[i]->loader.SetBodyAppliedForce(-fluid_density[i]/rho*system->GetBodyAppliedForce()); */
      }
      else {
        /* elems_loads_volumetric[i]->loader.SetBodyAppliedForce((1-fluid_density[i]/rho)*system->GetBodyAppliedForce()); */
      }
    }
  }
};

void cppCable::addNodestoContactCloud(std::shared_ptr<ChContactSurfaceNodeCloud> cloud) {
  for (int i = 0; i < nb_nodes-1; ++i) {
    if (beam_type == "BeamEuler") {
      cloud->AddNode(nodesRot[i], d);
    }
    else if (beam_type == "CableANCF") {
      cloud->AddNode(nodes[i], d);
    }
  }
};

cppMultiSegmentedCable * newMoorings(std::shared_ptr<ChSystem> system,
                                     std::shared_ptr<ChMesh> mesh,
                                     std::vector<double> length,
                                     std::vector<int> nb_elems,
                                     std::vector<double> d,
                                     std::vector<double> rho,
                                     std::vector<double> E,
                                     std::string beam_type)
{
  return new cppMultiSegmentedCable(system,
                                    mesh,
                                    length,
                                    nb_elems,
                                    d,
                                    rho,
                                    E,
                                    beam_type);
}


void cppAttachNodeToNodeFEAxyzD(cppMultiSegmentedCable* cable1,
                                int node1,
                                cppMultiSegmentedCable* cable2,
                                int node2) {
  auto con1 = chrono_types::make_shared<ChLinkNodeNode>();
  auto nodeA = cable1->nodes[node1];
  auto nodeB = cable2->nodes[node2];
  con1->Initialize(nodeA, nodeB);
  cable1->system->Add(con1);
}

void cppAttachNodeToNodeFEAxyzrot(cppMultiSegmentedCable* cable1,
                                  int node1,
                                  cppMultiSegmentedCable* cable2,
                                  int node2) {
  auto con1 = chrono_types::make_shared<ChLinkMateSpherical>();
  auto nodeA = cable1->nodesRot[node1];
  auto nodeB = cable2->nodesRot[node2];
  con1->Initialize(nodeA, nodeB, false, nodeA->GetPos(), nodeA->GetPos());
  cable1->system->Add(con1);
}
