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

  void ChBodyAddedMass::SetInertiaXX(const ChVector<>& iner) {
    variables.GetBodyInertia()(0, 0) = iner.x();
    variables.GetBodyInertia()(1, 1) = iner.y();
    variables.GetBodyInertia()(2, 2) = iner.z();
    variables.GetBodyInvInertia() = variables.GetBodyInertia().inverse();
    ChBody::SetInertiaXX(iner);
  }

  void ChBodyAddedMass::SetInertiaXY(const ChVector<>& iner) {
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
        assert(Mfullmass_in.GetRows() == variables.Get_ndof());
        assert(Mfullmass_in.GetColumns() == variables.Get_ndof());
        variables.SetMfullmass(Mfullmass_in);
    }

    // void ChBodyAddedMass::SetInvMfullmass(ChMatrixDynamic<> inv_Mfullmass_in) {
    //     assert(inv_Mfullmass_in.GetRows() == variables.Get_ndof());
    //     assert(inv_Mfullmass_in.GetColumns() == variables.Get_ndof());
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
    this->variables.Get_qb() = v.segment(off_v, 6);
    this->variables.Get_fb() = R.segment(off_v, 6);
}

void ChBodyAddedMass::IntFromDescriptor(const unsigned int off_v,  // offset in v
                               ChStateDelta& v,
                               const unsigned int off_L,  // offset in L
                               ChVectorDynamic<>& L) {
    v.segment(off_v, 6) = this->variables.Get_qb();
}

////

void ChBodyAddedMass::InjectVariables(ChSystemDescriptor& mdescriptor) {
    this->variables.SetDisabled(!this->IsActive());
    mdescriptor.InsertVariables(&this->variables);
}

void ChBodyAddedMass::VariablesFbReset() {
    this->variables.Get_fb().setZero();
}

void ChBodyAddedMass::VariablesFbLoadForces(double factor) {
    // add applied forces to 'fb' vector
    this->variables.Get_fb().segment(0, 3) += factor * Xforce.eigen();

    // add applied torques to 'fb' vector, including gyroscopic torque
    if (this->GetNoGyroTorque())
        this->variables.Get_fb().segment(3, 3) += factor * Xtorque.eigen();
    else
        this->variables.Get_fb().segment(3, 3) += factor * (Xtorque - gyro).eigen();
}

void ChBodyAddedMass::VariablesFbIncrementMq() {
    this->variables.Compute_inc_Mb_v(this->variables.Get_fb(), this->variables.Get_qb());
}

void ChBodyAddedMass::VariablesQbLoadSpeed() {
    // set current speed in 'qb', it can be used by the solver when working in incremental mode
    this->variables.Get_qb().segment(0, 3) = GetCoord_dt().pos.eigen();
    this->variables.Get_qb().segment(3, 3) = GetWvel_loc().eigen();
}

void ChBodyAddedMass::VariablesQbSetSpeed(double step) {
    ChCoordsys<> old_coord_dt = this->GetCoord_dt();

    // from 'qb' vector, sets body speed, and updates auxiliary data
    this->SetPos_dt(this->variables.Get_qb().segment(0, 3));
    this->SetWvel_loc(this->variables.Get_qb().segment(3, 3));

    // apply limits (if in speed clamping mode) to speeds.
    ClampSpeed();

    // compute auxiliary gyroscopic forces
    ComputeGyro();

    // Compute accel. by BDF (approximate by differentiation);
    if (step) {
        this->SetPos_dtdt((this->GetCoord_dt().pos - old_coord_dt.pos) / step);
        this->SetRot_dtdt((this->GetCoord_dt().rot - old_coord_dt.rot) / step);
    }
}

void ChBodyAddedMass::VariablesQbIncrementPosition(double dt_step) {
    if (!this->IsActive())
        return;

    // Updates position with incremental action of speed contained in the
    // 'qb' vector:  pos' = pos + dt * speed   , like in an Eulero step.

    ChVector<> newspeed(variables.Get_qb().segment(0, 3));
    ChVector<> newwel(variables.Get_qb().segment(3, 3));

    // ADVANCE POSITION: pos' = pos + dt * vel
    this->SetPos(this->GetPos() + newspeed * dt_step);

    // ADVANCE ROTATION: rot' = [dt*wwel]%rot  (use quaternion for delta rotation)
    ChQuaternion<> mdeltarot;
    ChQuaternion<> moldrot = this->GetRot();
    ChVector<> newwel_abs = Amatrix * newwel;
    double mangle = newwel_abs.Length() * dt_step;
    newwel_abs.Normalize();
    mdeltarot.Q_from_AngAxis(mangle, newwel_abs);
    ChQuaternion<> mnewrot = mdeltarot % moldrot;
    this->SetRot(mnewrot);
}

void ChBodyAddedMass::IntLoadResidual_F(const unsigned int off,  // offset in R residual
                               ChVectorDynamic<>& R,    // result: the R residual, R += c*F
                               const double c           // a scaling factor
                               ) {

    // add applied forces to 'fb' vector
    R.segment(off, 3) += c * Xforce.eigen();

    // add applied torques to 'fb' vector, including gyroscopic torque
    if (this->GetNoGyroTorque())
        R.segment(off + 3, 3) += c * Xtorque.eigen();
    else
        R.segment(off + 3, 3) += c * (Xtorque - gyro).eigen();
}

void ChBodyAddedMass::IntLoadResidual_Mv(const unsigned int off,      // offset in R residual
                                ChVectorDynamic<>& R,        // result: the R residual, R += c*M*v
                                const ChVectorDynamic<>& w,  // the w vector
                                const double c               // a scaling factor
                                ) {
  R(off + 0) += c * GetMass() * w(off + 0);
  R(off + 1) += c * GetMass() * w(off + 1);
  R(off + 2) += c * GetMass() * w(off + 2);
  ChVector<> Iw = GetInertia() * ChVector<>(w.segment(off + 3, 3));
  Iw *= c;
  R.segment(off + 3, 3) += Iw.eigen();
  ChMatrixDynamic<> ww = ChMatrixDynamic<>(6, 1);
  for (int i=0; i < 6; i++) {
    ww(i, 0) = w(off+i);
  }
  R.segment(off, 6) = variables.GetMfullmass()*ww*c;
}
}  // end namespace chrono

chrono::ChBodyAddedMass * newChBodyAddedMass()
{
  return new chrono::ChBodyAddedMass();
};
