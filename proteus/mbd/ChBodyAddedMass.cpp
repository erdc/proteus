#include "ChBodyAddedMass.h"
#include "ChVariablesBodyAddedMass.h"

namespace chrono {

  // Register into the object factory, to enable run-time dynamic creation and persistence
  CH_FACTORY_REGISTER(ChBodyAddedMass)

  ChBodyAddedMass::ChBodyAddedMass() {
    ChBody::variables = variables;
  }

  void ChBodyAddedMass::SetMass(double newmass) {
    variables.SetBodyMass(newmass);
    ChBody::variables.SetBodyMass(newmass);
  }

  void ChBodyAddedMass::SetInertia(const ChMatrix33<>& newXInertia) {
    variables.SetBodyInertia(newXInertia);
    ChBody::variables.SetBodyInertia(newXInertia);
  }

  void ChBodyAddedMass::SetInertiaXX(const ChVector3d& iner) {
    variables.GetBodyInertia()(0, 0) = iner.x();
    variables.GetBodyInertia()(1, 1) = iner.y();
    variables.GetBodyInertia()(2, 2) = iner.z();
    variables.GetBodyInvInertia() = variables.GetBodyInertia().inverse();
    ChBody::SetInertiaXX(iner);
  }

  void ChBodyAddedMass::SetInertiaXY(const ChVector3d& iner) {
    variables.GetBodyInertia()(0, 1) = iner.x();
    variables.GetBodyInertia()(0, 2) = iner.y();
    variables.GetBodyInertia()(1, 2) = iner.z();
    variables.GetBodyInertia()(1, 0) = iner.x();
    variables.GetBodyInertia()(2, 0) = iner.y();
    variables.GetBodyInertia()(2, 1) = iner.z();
    variables.GetBodyInvInertia() = variables.GetBodyInertia().inverse();
    ChBody::SetInertiaXY(iner);
  }

  //   ChVector<> ChBodyAddedMass::GetInertiaXX() {
  //       ChVector<> iner;
  //       iner.x() = variables.GetBodyInertia().GetElement(0, 0);
  //       iner.y() = variables.GetBodyInertia().GetElement(1, 1);
  //       iner.z() = variables.GetBodyInertia().GetElement(2, 2);
  //       return iner;
  //   }

  //   ChVector<> ChBodyAddedMass::GetInertiaXY() {
  //       ChVector<> iner;
  //       iner.x() = variables.GetBodyInertia().GetElement(0, 1);
  //       iner.y() = variables.GetBodyInertia().GetElement(0, 2);
  //       iner.z() = variables.GetBodyInertia().GetElement(1, 2);
  //       return iner;
  //   }

  void ChBodyAddedMass::SetMfullmass(ChMatrixDynamic<> Mfullmass_in) {
    assert(Mfullmass_in.rows() == variables.GetDOF());
    assert(Mfullmass_in.cols() == variables.GetDOF());
    variables.SetMfullmass(Mfullmass_in);
  }

  // void ChBodyAddedMass::SetInvMfullmass(ChMatrixDynamic<> inv_Mfullmass_in) {
  //     assert(inv_Mfullmass_in.GetRows() == variables.GetDOF());
  //     assert(inv_Mfullmass_in.GetColumns() == variables.GetDOF());
  //     ChMatrixDynamic<>& Mm = variables.GetInvMfullmass();
  //     for (int i = 0; i < 6; i++) {
  //         for (int j = 0; j < 6; j++) {
  //             Mm(i, j) = inv_Mfullmass_in(i, j);
  //         }
  //     }
  // }
  //// STATE BOOKKEEPING FUNCTIONS



  void ChBodyAddedMass::IntToDescriptor(const unsigned int off_v,  // offset in v, R
					const ChStateDelta& v,
					const ChVectorDynamic<>& R,
					const unsigned int off_L,  // offset in L, Qc
					const ChVectorDynamic<>& L,
					const ChVectorDynamic<>& Qc) {
    this->variables.State() = v.segment(off_v, 6);
    this->variables.Force() = R.segment(off_v, 6);
  }

  void ChBodyAddedMass::IntFromDescriptor(const unsigned int off_v,  // offset in v
					  ChStateDelta& v,
					  const unsigned int off_L,  // offset in L
					  ChVectorDynamic<>& L) {
    v.segment(off_v, 6) = this->variables.State();
  }

  ////

  void ChBodyAddedMass::InjectVariables(ChSystemDescriptor& mdescriptor) {
    this->variables.SetDisabled(!this->IsActive());
    mdescriptor.InsertVariables(&this->variables);
  }

  void ChBodyAddedMass::VariablesFbReset() {
    this->variables.Force().setZero();
  }

  void ChBodyAddedMass::VariablesFbLoadForces(double factor) {
    // add applied forces to Force vector
    this->variables.Force().segment(0, 3) += factor * Xforce.eigen();

    // add applied torques to Force vector, including gyroscopic torque
    if (this->IsUsingGyroTorque())
      this->variables.Force().segment(3, 3) += factor * Xtorque.eigen();
    else
      this->variables.Force().segment(3, 3) += factor * (Xtorque - gyro).eigen();
  }

  void ChBodyAddedMass::VariablesFbIncrementMq() {
    this->variables.AddMassTimesVector(this->variables.Force(), this->variables.State());
  }

  void ChBodyAddedMass::VariablesQbLoadSpeed() {
    // set current speed in State, it can be used by the solver when working in incremental mode
    this->variables.State().segment(0, 3) = GetCoordsysDt().pos.eigen();
    this->variables.State().segment(3, 3) = GetAngVelLocal().eigen();
  }

  void ChBodyAddedMass::VariablesQbSetSpeed(double step) {
    ChCoordsys<> old_coord_dt = this->GetCoordsysDt();

    // from State vector, sets body speed, and updates auxiliary data
    this->SetPosDt(this->variables.State().segment(0, 3));
    this->SetAngVelLocal(this->variables.State().segment(3, 3));

    // apply limits (if in speed clamping mode) to speeds.
    ClampSpeed();

    // compute auxiliary gyroscopic forces
    ComputeGyro();

    // Compute accel. by BDF (approximate by differentiation);
    if (step) {
      this->SetPosDt2((this->GetCoordsysDt().pos - old_coord_dt.pos) / step);
      this->SetRotDt2((this->GetCoordsysDt().rot - old_coord_dt.rot) / step);
    }
  }

  void ChBodyAddedMass::VariablesQbIncrementPosition(double dt_step) {
    if (!this->IsActive())
      return;

    // Updates position with incremental action of speed contained in the
    // State vector:  pos' = pos + dt * speed   , like in an Eulero step.

    ChVector3d newspeed(variables.State().segment(0, 3));
    ChVector3d newwel(variables.State().segment(3, 3));

    // ADVANCE POSITION: pos' = pos + dt * vel
    this->SetPos(this->GetPos() + newspeed * dt_step);

    // ADVANCE ROTATION: rot' = [dt*wwel]%rot  (use quaternion for delta rotation)
    ChQuaternion<> mdeltarot;
    ChQuaternion<> moldrot = this->GetRot();
    ChVector3d newwel_abs = this->GetRotMat() * newwel;
    double mangle = newwel_abs.Length() * dt_step;
    newwel_abs.Normalize();
    mdeltarot.SetFromAngleAxis(mangle, newwel_abs);
    ChQuaternion<> mnewrot = mdeltarot * moldrot;
    this->SetRot(mnewrot);
  }

  void ChBodyAddedMass::IntLoadResidual_F(const unsigned int off,  // offset in R residual
					  ChVectorDynamic<>& R,    // result: the R residual, R += c*F
					  const double c           // a scaling factor
					  ) {

    // add applied forces to Force vector
    R.segment(off, 3) += c * Xforce.eigen();

    // add applied torques to Force vector, including gyroscopic torque
    if (this->IsUsingGyroTorque())
      R.segment(off + 3, 3) += c * Xtorque.eigen();
    else
      R.segment(off + 3, 3) += c * (Xtorque - gyro).eigen();
  }

  void ChBodyAddedMass::IntLoadResidual_Mv(const unsigned int off,      // offset in R residual
					   ChVectorDynamic<>& R,        // result: the R residual, R += c*M*v
					   const ChVectorDynamic<>& w,  // the w vector
					   const double c               // a scaling factor
					   ) {
    ChMatrixDynamic<> ww = ChMatrixDynamic<>(6, 1);
    for (int i=0; i < 6; i++) {
      ww(i, 0) = w(off+i);
    }
    R.segment(off, 6) += variables.GetMfullmass()*ww*c;
  }
}  // end namespace chrono

chrono::ChBodyAddedMass * newChBodyAddedMass()
{
  return new chrono::ChBodyAddedMass();
};
