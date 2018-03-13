//#pragma once

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/core/ChTransform.h"
#include <iostream>
#include <fstream>

using namespace chrono;
using namespace std;



class cppSystem {
 public:
  ChSystemSMC system;
  double* gravity;
  double chrono_dt;
  std::string directory;
  cppSystem(double* gravity);
  void step(double proteus_dt, int n_substeps);
  void setChTimeStep(double dt);
  void setGravity(double* gravity);
  void setDirectory(std::string dir);
  void setTimestepperType(std::string tstype, bool verbose);
  void setCollisionEnvelopeMargin(double envelope, double margin);
};


class cppRigidBody {
 public:
  ChVector<> free_x;
  ChVector<> free_r;
  ChVector<> pos;
  ChVector<> pos_last;
  ChVector<> vel;
  ChVector<> vel_last;
  ChVector<> acc;
  ChVector<> acc_last;
  ChVector<> angvel;
  ChVector<> angvel_last;
  ChVector<> angacc;
  ChVector<> angacc_last;
  ChMatrix33<double> rotm;
  ChMatrix33<double> rotm_last;
  ChQuaternion<double> rotq;
  ChQuaternion<double> rotq_last;
  ChVector<> F;
  ChVector<> F_last;
  ChVector<> M;
  ChVector<> M_last;
  std::shared_ptr<ChLinkLockLock> lock_motion;  // lock for prescribed motion
  double lock_motion_t_max;  // max time up to which lock motion is used
  double mass;
  double mooring_restlength;
  std::shared_ptr<ChLinkSpring> spring;
  /* ChVector <> inertia; */
  double* inertia;
  std::shared_ptr<ChBody> body;
  cppSystem* system;
  cppRigidBody(cppSystem* system);
  ChVector<double> hxyz(double* x, double t);
  double hx(double* x, double t);
  double hy(double* x, double t);
  double hz(double* x, double t);
  void prestep(double* force, double* torque);
  void poststep();
  void setRotation(double* quat);
  void setPosition(double* quat);
  void setConstraints(double* free_x, double* free_y);
  void setInertiaXX(double* inertia);
  void addSpring(double stiffness,
                 double damping,
                 double* fairlead,
                 double* anchor,
                 double rest_length);
  void addPrismaticLinksWithSpring(double* pris1,
                                   double* pris2,
                                   double stiffness,
                                   double damping,
                                   double rest_length);
  void addPrismaticLinkX(double* pris1);
  void setName(std::string name);
  void setPrescribedMotionPoly(double coeff1);
  void setPrescribedMotionSine(double a, double f);
  void setPrescribedMotionCustom(std::vector<double> t, std::vector<double> x,
                                 std::vector<double> y, std::vector<double> z,
                                 std::vector<double> ang, std::vector<double> ang2,
                                 std::vector<double> ang3, double t_max);
};

cppSystem::cppSystem(double* gravity):
gravity(gravity)
{
  chrono_dt = 0.000001;
  system.Set_G_acc(ChVector<>(gravity[0], gravity[1], gravity[2]));
  directory = "./";
  // SOLVER OPTIONS
  system.SetSolverType(ChSolver::Type::MINRES);  // SOLVER_MINRES: good convergence, supports FEA, does not support DVI yet
  auto msolver = std::static_pointer_cast<ChSolverMINRES>(system.GetSolver());
  msolver->SetDiagonalPreconditioning(true);
  system.SetSolverWarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
  system.SetMaxItersSolverSpeed(100); // max iteration for iterative solvers
  system.SetMaxItersSolverStab(100); // max iteration for stabilization (iterative solvers)
  system.SetTolForce(1e-10);
  //system.SetMaxItersSolverSpeed(100);  
  //system.SetMaxItersSolverStab(100);  
  //system.SetTolForce(1e-14); // default: 0.001
  //system.SetMaxiter(200); // default: 6. Max constraints to reach tolerance on constraints.
  //system.SetTol(1e-10); // default: 0.0002. Tolerance for keeping constraints together.
  system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED); // used before: ChSystemSMC::INT_EULER_IMPLICIT_LINEARIZED
  if (auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(system.GetTimestepper())) {
    mystepper->SetAlpha(-0.2);
  }
}

void cppSystem::setTimestepperType(std::string tstype, bool verbose=false) {
  if (tstype == "HHT") {
    system.SetTimestepperType(ChTimestepper::Type::HHT);
      auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(system.GetTimestepper());
      mystepper->SetAlpha(-0.2);
      mystepper->SetMaxiters(10);
      mystepper->SetAbsTolerances(1e-6);
      mystepper->SetMode(ChTimestepperHHT::POSITION);
      mystepper->SetScaling(false);
      mystepper->SetVerbose(verbose);
      mystepper->SetModifiedNewton(false);
    }
    else if (tstype == "Euler") {
      system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
    }
    else if (tstype == "Trapezoidal") {
      system.SetTimestepperType(ChTimestepper::Type::TRAPEZOIDAL);
    }
  }

void cppSystem::setGravity(double* gravity)
{
  system.Set_G_acc(ChVector<>(gravity[0], gravity[1], gravity[2]));
}

void cppSystem::step(double proteus_dt, int n_substeps=1)
{
    double dt2 = proteus_dt/(double)n_substeps;
   for (int i = 0; i < n_substeps; ++i) {
     system.DoStepDynamics(dt2);
   }
}

void cppSystem::setChTimeStep(double dt) {
    chrono_dt = dt;
};

cppRigidBody::cppRigidBody(cppSystem* system):
  system(system)
{

  body = std::make_shared<ChBody>();
  // add body to system
  /* system->system.AddBody(body); */ // now added externally in cython
  // basic attributes of body
  rotm = body->GetA();
  rotm_last = body->GetA();
  pos = body->GetPos();
  pos_last = body->GetPos();
  body->SetMass(mass);
  free_x = {1., 1., 1.};
  free_r = {1., 1., 1.};
  lock_motion_t_max = 0.;
}


void cppSystem::setDirectory(std::string dir) {
    directory = dir;
}

ChVector<double> cppRigidBody::hxyz(double* x, double t)
{
  /* rotm = body->GetA(); */
  ChVector<double> xx = ChVector<double>(x[0], x[1], x[2]);
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(xx, pos_last, rotq_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotq);
  return xNew - xx;
}


void cppSystem::setCollisionEnvelopeMargin(double envelope, double margin) {
  collision::ChCollisionModel::SetDefaultSuggestedEnvelope(envelope);
  collision::ChCollisionModel::SetDefaultSuggestedMargin(margin);
}

double cppRigidBody::hx(double* x, double t)
{
  /* rotm = body->GetA(); */
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotq_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotq);
  return xNew.x() - x[0];
}

double cppRigidBody::hy(double* x, double t)
{
  /* rotm = body->GetA(); */
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotq_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotq);
  return xNew.y() - x[1];
}

double cppRigidBody::hz(double* x, double t)
{
  /* rotm = body->GetA(); */
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotq_last);
  ChVector<double> xNew = ChTransform<double>::TransformLocalToParent(local, pos, rotq);
  return xNew.z() - x[2];
}

void cppRigidBody::prestep(double* force, double* torque)
{
  /* step to call before running chrono system step */
  pos_last = body->GetPos();
  vel_last = body->GetPos_dt();
  acc_last = body->GetPos_dtdt();
  rotm_last = body->GetA();
  rotq_last = body->GetRot();
  angacc_last = body->GetWacc_loc();
  angvel_last = body->GetWvel_loc();
  F_last = body->Get_Xforce();
  M_last = body->Get_Xtorque();
  // apply external forces
  body->Empty_forces_accumulators();
  // calculate opposite force of gravity if free_x is 0
  double forceG[3]={0.,0.,0.};
  if (free_x.x() == 0) {forceG[0] = -system->system.Get_G_acc().x()*body->GetMass();}
  if (free_x.y() == 0) {forceG[1] = -system->system.Get_G_acc().y()*body->GetMass();}
  if (free_x.z() == 0) {forceG[2] = -system->system.Get_G_acc().z()*body->GetMass();}
  body->Accumulate_force(ChVector<double>(forceG[0]+force[0]*free_x.x(),
                                          forceG[1]+force[1]*free_x.y(),
                                          forceG[2]+force[2]*free_x.z()),
                         pos_last,
                         false);
  body->Accumulate_torque(ChVector<double>(torque[0]*free_r.x(),
                                           torque[1]*free_r.y(),
                                           torque[2]*free_r.z()),
                          true);
  if (spring!=0) {
      double spring_length = spring->Get_SpringLength();
      if (spring_length < mooring_restlength) {
          spring->SetDisabled(true);//Set_SpringRestLength(spring_length);
      }
      else {
          spring->SetDisabled(false);//Set_SpringRestLength(mooring_restlength);
      }
  }
}



void cppRigidBody::poststep()
{
  pos = body->GetPos();
  vel = body->GetPos_dt();
  acc = body->GetPos_dtdt();
  rotm = body->GetA();
  rotq = body->GetRot();
  angacc = body->GetWacc_loc();
  angvel = body->GetWvel_loc();
  F = body->Get_Xforce();
  M = body->Get_Xtorque();
  if (lock_motion_t_max > 0) {
    double t = system->system.GetChTime();
    if (lock_motion_t_max < t && lock_motion->IsDisabled() == false) {
      lock_motion->SetDisabled(true);
    }
  }
}

void cppRigidBody::setPrescribedMotionCustom(std::vector<double> t,
                                             std::vector<double> x,
                                             std::vector<double> y,
                                             std::vector<double> z,
                                             std::vector<double> ang,
                                             std::vector<double> ang2,
                                             std::vector<double> ang3,
                                             double t_max) {
  auto fixed_body = std::make_shared<ChBody>();
  fixed_body->SetPos(body->GetPos());
  fixed_body->SetBodyFixed(true);
  system->system.Add(fixed_body);
  lock_motion = std::make_shared<ChLinkLockLock>();
  lock_motion_t_max = t_max;
  lock_motion->Initialize(body, fixed_body, fixed_body->GetCoord());
  system->system.Add(lock_motion);
  if (x.size() > 0) {
    auto forced_motion = std::make_shared<ChFunction_Recorder>();
    for (int i = 0; i < x.size(); i++) {
      forced_motion->AddPoint(t[i], x[i]);
    }
    std::shared_ptr<ChFunction> forced_ptr = forced_motion;
    lock_motion->SetMotion_X(forced_ptr);
  }
  if (y.size() > 0) {
    auto forced_motion = std::make_shared<ChFunction_Recorder>();
    for (int i = 0; i < y.size(); i++) {
      forced_motion->AddPoint(t[i], y[i]);
    }
    std::shared_ptr<ChFunction> forced_ptr = forced_motion;
    lock_motion->SetMotion_Y(forced_ptr);
  }
  if (z.size() > 0) {
    auto forced_motion = std::make_shared<ChFunction_Recorder>();
    for (int i = 0; i < z.size(); i++) {
      forced_motion->AddPoint(t[i], z[i]);
    }
    std::shared_ptr<ChFunction> forced_ptr = forced_motion;
    lock_motion->SetMotion_Z(forced_ptr);
  }
  if (ang.size() > 0) {
    auto forced_motion = std::make_shared<ChFunction_Recorder>();
    for (int i = 0; i < ang.size(); i++) {
      forced_motion->AddPoint(t[i], ang[i]);
    }
    std::shared_ptr<ChFunction> forced_ptr = forced_motion;
    lock_motion->SetMotion_ang(forced_ptr);
  }
  if (ang2.size() > 0) {
    auto forced_motion = std::make_shared<ChFunction_Recorder>();
    for (int i = 0; i < ang2.size(); i++) {
      forced_motion->AddPoint(t[i], ang2[i]);
    }
    std::shared_ptr<ChFunction> forced_ptr = forced_motion;
    lock_motion->SetMotion_ang2(forced_ptr);
  }
  if (ang3.size() > 0) {
    auto forced_motion = std::make_shared<ChFunction_Recorder>();
    for (int i = 0; i < ang3.size(); i++) {
      forced_motion->AddPoint(t[i], ang3[i]);
    }
    std::shared_ptr<ChFunction> forced_ptr = forced_motion;
    lock_motion->SetMotion_ang3(forced_ptr);
  }
}

void cppRigidBody::setPrescribedMotionPoly(double coeff1) {
  auto fixed_body = std::make_shared<ChBody>();
  fixed_body->SetPos(body->GetPos());
  fixed_body->SetBodyFixed(true);
  system->system.Add(fixed_body);
  auto lock = std::make_shared<ChLinkLockLock>();
  lock->Initialize(body, fixed_body, fixed_body->GetCoord());
  system->system.Add(lock);
  auto forced_motion = std::make_shared<ChFunction_Poly>();
  forced_motion->Set_order(1);
  forced_motion->Set_coeff(coeff1, 1);
  std::shared_ptr<ChFunction> forced_ptr = forced_motion;
  lock->SetMotion_X(forced_ptr);
}


void cppRigidBody::setPrescribedMotionSine(double a, double f) {
  auto fixed_body = std::make_shared<ChBody>();
  fixed_body->SetPos(body->GetPos());
  fixed_body->SetBodyFixed(true);
  system->system.Add(fixed_body);
  auto lock = std::make_shared<ChLinkLockLock>();
  lock->Initialize(body, fixed_body, fixed_body->GetCoord());
  system->system.Add(lock);
  auto forced_motion = std::make_shared<ChFunction_Sine>();
  forced_motion->Set_amp(a);
  forced_motion->Set_freq(f);
  std::shared_ptr<ChFunction> forced_ptr = forced_motion;
  lock->SetMotion_X(forced_ptr);
}

void cppRigidBody::setPosition(double* position){
  body->SetPos(ChVector<>(position[0], position[1], position[2]));
}

void cppRigidBody::setRotation(double* quat) {
  body->SetRot(ChQuaternion<double>(quat[0], quat[1], quat[2], quat[3]));
}

void cppRigidBody::setConstraints(double* free_x_in, double* free_r_in){
  free_x = ChVector<>(free_x_in[0], free_x_in[1], free_x_in[2]);
  free_r = ChVector<>(free_r_in[0], free_r_in[1], free_r_in[2]);
}

void cppRigidBody::setInertiaXX(double* inertia){
  body->SetInertiaXX(ChVector<>(inertia[0], inertia[1], inertia[2]));
}


void cppRigidBody::addSpring(double stiffness,
                             double damping,
                             double* fairlead,
                             double* anchor,
                             double rest_length)
{
  mooring_restlength = rest_length;
  spring = std::make_shared<ChLinkSpring>();
  std::shared_ptr<ChBody> anchor_body = std::make_shared<ChBody>();
  anchor_body->SetPos(ChVector<>(anchor[0], anchor[1], anchor[2]));
  anchor_body->SetBodyFixed(true);
  system->system.AddBody(anchor_body);
  spring->Initialize(body,
                     anchor_body,
                     true, // true for pos relative to bodies
                     ChVector<>(fairlead[0], fairlead[1], fairlead[2]),
                     ChVector<>(0.,0.,0.),
                     false,  // true for auto rest length (distance between body1 and body2)
                     rest_length);
  spring->Set_SpringK(stiffness);
  spring->Set_SpringR(damping);
  system->system.AddLink(spring);
}

void cppRigidBody::addPrismaticLinkX(double* pris1)
{
  auto mybod2 = std::make_shared<ChBody>();
  mybod2->SetName("PRIS1");
  mybod2->SetPos(ChVector<>(pris1[0], pris1[1], pris1[2]));
  mybod2->SetMass(0.00001);
  mybod2->SetBodyFixed(true);
  system->system.AddBody(mybod2);
  auto mylink1 = std::make_shared<ChLinkLockPrismatic>();
  auto mycoordsys1 = ChCoordsys<>(mybod2->GetPos(),Q_from_AngAxis(CH_C_PI/2., VECT_Y));//Q_from_AngAxis(CH_C_PI / 2, VECT_X));
  mylink1->Initialize(mybod2, body, mycoordsys1);
  system->system.AddLink(mylink1);
}

void cppRigidBody::addPrismaticLinksWithSpring(double* pris1,
                                               double* pris2,
                                               double stiffness,
                                               double damping,
                                               double rest_length)
{
  mooring_restlength = rest_length;
  auto fairlead = std::make_shared<ChBody>();
  fairlead->SetName("PRIS3");
  fairlead->SetPos(body->GetPos());
  fairlead->SetMass(0.00001);
  system->system.AddBody(fairlead);
  auto mybod2 = std::make_shared<ChBody>();
  mybod2->SetName("PRIS1");
  mybod2->SetPos(ChVector<>(pris1[0], pris1[1], pris1[2]));
  mybod2->SetMass(0.00001);
  //mybod2->AddForce(-system->system.Get_G_acc());
  //mybod2->SetBodyFixed(true);
  system->system.AddBody(mybod2);
  auto mybod3 = std::make_shared<ChBody>();
  mybod3->SetName("PRIS2");
  mybod3->SetPos(ChVector<>(pris2[0], pris2[1], pris2[2]));
  mybod3->SetBodyFixed(true);
  system->system.AddBody(mybod3);

  auto mylink1 = std::make_shared<ChLinkLockPrismatic>();
  system->system.AddLink(mylink1);
  auto mycoordsys1 = ChCoordsys<>(mybod2->GetPos(),Q_from_AngAxis(CH_C_PI/2., VECT_Y));//Q_from_AngAxis(CH_C_PI / 2, VECT_X));
  mylink1->Initialize(fairlead, mybod2, mycoordsys1);



  auto mylink2 = std::make_shared<ChLinkLockPrismatic>();
  system->system.AddLink(mylink2);
  auto mycoordsys2 = ChCoordsys<>(mybod3->GetPos(),Q_from_AngAxis(CH_C_PI/2., VECT_X));//Q_from_AngAxis(CH_C_PI / 2, VECT_X));
  mylink2->Initialize(mybod2, mybod3,mycoordsys2);

  auto mylink3 = std::make_shared<ChLinkLockSpherical>();
  //auto mylink3 = std::make_shared<ChLinkLockRevolute>();
  //mylink3->SetMotion_axis(ChVector<>(0.,1.,0.));
  system->system.AddLink(mylink3);
  mylink3->Initialize(fairlead, body, false, fairlead->GetCoord(), body->GetCoord());



  spring = std::make_shared<ChLinkSpring>();
  spring->Initialize(fairlead,
                     mybod2,
                     true, // true for pos relative to bodies
                     ChVector<>(0.,0.,0.),
                     ChVector<>(0.,0.,0.),
                     false,
                     rest_length);  // true for auto rest length (distance between body1 and body2));
  spring->Set_SpringK(stiffness);
  spring->Set_SpringR(damping);
  spring->SetName("SPRING1");
  system->system.AddLink(spring);
}

void cppRigidBody::setName(std::string name) {
  body->SetNameString(name);
}


cppSystem * newSystem(double* gravity)
{
  return new cppSystem(gravity);
}



cppRigidBody * newRigidBody(cppSystem* system)
{
  return new cppRigidBody(system);
}



void ChLinkLockBodies(std::shared_ptr<ChBody> body1,
                      std::shared_ptr<ChBody> body2,
                      ChSystemSMC& system,
                      ChCoordsys<> coordsys,
                      double limit_X=0.,
                      double limit_Y=0.,
                      double limit_Z=0.,
                      double limit_Rx=0.,
                      double limit_Ry=0.,
                      double limit_Rz=0.) {
  auto mylink = std::make_shared<ChLinkLock>();
  system.AddLink(mylink);
  auto chlimit_X = ChLinkLimit();
  chlimit_X.Set_active(true);
  chlimit_X.Set_max(limit_X);
  auto chlimit_Y = ChLinkLimit();
  chlimit_X.Set_active(true);
  chlimit_Y.Set_max(limit_Y);
  auto chlimit_Z = ChLinkLimit();
  chlimit_X.Set_active(true);
  chlimit_Z.Set_max(limit_Z);
  auto chlimit_Rx = ChLinkLimit();
  chlimit_Rx.Set_max(limit_Rx);
  chlimit_X.Set_active(true);
  auto chlimit_Ry = ChLinkLimit();
  chlimit_X.Set_active(true);
  chlimit_Ry.Set_max(limit_Ry);
  auto chlimit_Rz = ChLinkLimit();
  chlimit_X.Set_active(true);
  chlimit_Rz.Set_max(limit_Rz);
  mylink->SetLimit_X(&chlimit_X);
  mylink->SetLimit_Y(&chlimit_Y);
  mylink->SetLimit_Z(&chlimit_Z);
  mylink->SetLimit_Rx(&chlimit_Rx);
  mylink->SetLimit_Ry(&chlimit_Ry);
  mylink->SetLimit_Rz(&chlimit_Rz);
  mylink->Initialize(body1, body2, coordsys);
}
