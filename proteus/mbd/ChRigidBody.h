#include "chrono/physics/ChSystemDEM.h"
#include "chrono/timestepper/ChTimestepper.h"
#include "chrono/solver/ChSolverMINRES.h"
#include <memory>
#include <iostream>
#include <fstream>
using namespace chrono;
using namespace std;


class cppSystem {
 public:
  ChSystemDEM system;
  double* gravity;
  double chrono_dt;
  std::string directory;
  cppSystem(double* gravity);
  void step(double proteus_dt, int n_substeps);
  void setChTimeStep(double dt);
  void recordBodyList();
  void setGravity(double* gravity);
  void setDirectory(std::string dir);
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
  double mass;
  double mooring_restlength;
  std::shared_ptr<ChLinkSpring> spring;
  /* ChVector <> inertia; */
  double* inertia;
  std::shared_ptr<ChBody> body;
  cppSystem* system;
  cppRigidBody(cppSystem* system,
               double* pos,
               double* rotq,
               double mass,
               double* inertia,
               double* free_x,
               double* free_r);
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
  void setName(std::string name);
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
  system.SetMaxItersSolverSpeed(100);  // max iteration for iterative solvers
  system.SetMaxItersSolverStab(100);  // max iteration for stabilization (iterative solvers)
  system.SetTolForce(1e-14); // default: 0.001
  //system.SetMaxiter(200); // default: 6. Max constraints to reach tolerance on constraints.
  //system.SetTol(1e-10); // default: 0.0002. Tolerance for keeping constraints together.
  //system.SetIntegrationType(ChSystemDEM::INT_HHT); // used before: ChSystemDEM::INT_EULER_IMPLICIT_LINEARIZED
  system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED); // used before: ChSystemDEM::INT_EULER_IMPLICIT_LINEARIZED
  if (auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(system.GetTimestepper())) {
    mystepper->SetAlpha(-0.2);
  }
  //system.SetTimestepper(std::make_shared<ChTimestepperEulerImplicitLinearized>());  // default: fast, 1st order
}


void cppSystem::setGravity(double* gravity)
{
  system.Set_G_acc(ChVector<>(gravity[0], gravity[1], gravity[2]));
}

void cppSystem::step(double proteus_dt, int n_substeps=1)
{
  /* double substeps = proteus_dt/chrono_dt+0.5; //+0.5 for rounding n substeps */

  /* int n_substeps = (int)substeps; */
    double dt2 = proteus_dt/(double)n_substeps;
   for (int i = 0; i < n_substeps; ++i) {
     system.DoStepDynamics(dt2);
   }
  /* double t = chrono_dt; */
  /* if (t > proteus_dt) { */
  /*     system.DoStepDynamics(proteus_dt); */
  /*   } */
  /* else { */
  /*     while (t <= proteus_dt) { */
  /*         system.DoStepDynamics(chrono_dt); */
  /*         t += chrono_dt; */
  /*     } */
  /*     if (t != proteus_dt) {  // means t went above dt, need last time step */
  /*         system.DoStepDynamics(proteus_dt-(t-chrono_dt)); */
  /* } */
  /* } */
}

void cppSystem::recordBodyList() {
      std::vector<std::shared_ptr<ChBody>>& bodylist = *system.Get_bodylist();
      double t = system.GetChTime();
      if (t == 0) {
          for (int i = 0; i<bodylist.size(); i++) {
              std::shared_ptr<ChBody> bod = bodylist[i];
              fstream myfile;
              myfile.open (directory+bod->GetNameString()+".csv", std::ios_base::out);
              myfile << "t,x,y,z,e0,e1,e2,e3,ux,uy,uz,ax,ay,az,Fx,Fy,Fz,Mx,My,Mz,\n";
              myfile.close();
          }
      }
      for (int i = 0; i<bodylist.size(); i++) {
          std::shared_ptr<ChBody> bod = bodylist[i];
          fstream myfile;
          myfile.open (directory+bod->GetNameString()+".csv", std::ios_base::app);     
          ChVector<> bpos = bod->GetPos();
          ChVector<> bvel = bod->GetPos_dt();
          ChVector<> bacc = bod->GetPos_dtdt();
          ChVector<> bfor = bod->Get_Xforce();
          ChVector<> btor = bod->Get_Xtorque();
          ChQuaternion<> brot = bod->GetRot();
          myfile << t << ",";     
          myfile << bpos.x() << "," << bpos.y() << "," << bpos.z() << ",";     
          myfile << brot.e0() << "," << brot.e1() << "," << brot.e2() << "," << brot.e3() << ",";     
          myfile << bvel.x() << "," << bvel.y() << "," << bvel.z() << ",";     
          myfile << bacc.x() << "," << bacc.y() << "," << bacc.z() << ",";     
          myfile << bfor.x() << "," << bfor.y() << "," << bfor.z() << ",";     
          myfile << btor.x() << "," << btor.y() << "," << btor.z() << ",";     
          myfile << "\n";        
          myfile.close();
    }
}


void cppSystem::setChTimeStep(double dt) {
    chrono_dt = dt;
};

cppRigidBody::cppRigidBody(cppSystem* system,
                           double* posin,
                           double* rotin,
                           double mass,
                           double* inertia,
                           double* free_xin,
                           double* free_rin):
  system(system),
  mass(mass),
  inertia(inertia),
  free_x(free_xin[0], free_xin[1], free_xin[2]),
  free_r(free_rin[0], free_rin[1], free_rin[2])
{

  body = std::make_shared<ChBody>();
  // add body to system
  system->system.AddBody(body);
  // basic attributes of body
  pos = ChVector<>(posin[0], posin[1], posin[2]);
  rotq = ChQuaternion<>(rotin[0], rotin[1], rotin[2], rotin[3]);
  body->SetPos(pos);
  body->SetRot(rotq);
  body->SetInertiaXX(ChVector<>(1.,
                                1.,
                                inertia[2]));  // careful division by zero!
  if (free_x.x() == 0 && free_x.y() == 0 && free_x.z() ==0 && free_r.x() == 0 && free_r.y() ==0 && free_r.z() == 0) {
  body->SetBodyFixed(true);
  }
  rotm = body->GetA();
  rotm_last = body->GetA();
  pos = body->GetPos();
  pos_last = body->GetPos();
  body->SetMass(mass);
}


void cppSystem::setDirectory(std::string dir) {
    directory = dir;
}


double cppRigidBody::hx(double* x, double t)
{
  rotm = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotm_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotm);
  return xNew.x() - x[0];
}

double cppRigidBody::hy(double* x, double t)
{
  rotm = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotm_last);
  ChVector<double> xNew  = ChTransform<double>::TransformLocalToParent(local, pos, rotm);
  return xNew.y() - x[1];
}

double cppRigidBody::hz(double* x, double t)
{
  rotm = body->GetA();
  ChVector<double> local = ChTransform<double>::TransformParentToLocal(ChVector<double>(x[0],x[1],x[2]), pos_last, rotm_last);
  ChVector<double> xNew = ChTransform<double>::TransformLocalToParent(local, pos, rotm);
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
  angacc_last = body->GetWvel_loc();
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
  angacc = body->GetWvel_loc();
  angvel = body->GetWvel_loc();
  F = body->Get_Xforce();
  M = body->Get_Xtorque();
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



cppRigidBody * newRigidBody(cppSystem* system,
                            double* position,
                            double* rotq,
                            double mass,
                            double* inertia,
                            double* free_x,
                            double* free_r)
{
  return new cppRigidBody(system,
                          position,
                          rotq,
                          mass,
                          inertia,
                          free_x,
                          free_r);
}
