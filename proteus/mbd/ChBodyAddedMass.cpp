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
    variables.GetBodyInertia().SetElement(0, 0, iner.x());
    variables.GetBodyInertia().SetElement(1, 1, iner.y());
    variables.GetBodyInertia().SetElement(2, 2, iner.z());
    variables.GetBodyInertia().FastInvert(variables.GetBodyInvInertia());
    ChBody::SetInertiaXX(iner);
  }

  void ChBodyAddedMass::SetInertiaXY(const ChVector<>& iner) {
    variables.GetBodyInertia().SetElement(0, 1, iner.x());
    variables.GetBodyInertia().SetElement(0, 2, iner.y());
    variables.GetBodyInertia().SetElement(1, 2, iner.z());
    variables.GetBodyInertia().SetElement(1, 0, iner.x());
    variables.GetBodyInertia().SetElement(2, 0, iner.y());
    variables.GetBodyInertia().SetElement(2, 1, iner.z());
    variables.GetBodyInertia().FastInvert(variables.GetBodyInvInertia());
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
        ChMatrix<>& Mm = variables.GetMfullmass();
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                Mm.SetElement(i, j, Mfullmass_in.GetElement(i, j));
            }
        }
    }

    void ChBodyAddedMass::SetInvMfullmass(ChMatrixDynamic<> inv_Mfullmass_in) {
        assert(inv_Mfullmass_in.GetRows() == variables.Get_ndof());
        assert(inv_Mfullmass_in.GetColumns() == variables.Get_ndof());
        ChMatrix<>& Mm = variables.GetInvMfullmass();
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                Mm(i, j) = inv_Mfullmass_in.GetElement(i, j);
            }
        }
    }
//// STATE BOOKKEEPING FUNCTIONS



void ChBodyAddedMass::IntToDescriptor(const unsigned int off_v,  // offset in v, R
                             const ChStateDelta& v,
                             const ChVectorDynamic<>& R,
                             const unsigned int off_L,  // offset in L, Qc
                             const ChVectorDynamic<>& L,
                             const ChVectorDynamic<>& Qc) {
    this->variables.Get_qb().PasteClippedMatrix(v, off_v, 0, 6, 1, 0, 0);  // for solver warm starting only
    this->variables.Get_fb().PasteClippedMatrix(R, off_v, 0, 6, 1, 0, 0);  // solver known term
}

void ChBodyAddedMass::IntFromDescriptor(const unsigned int off_v,  // offset in v
                               ChStateDelta& v,
                               const unsigned int off_L,  // offset in L
                               ChVectorDynamic<>& L) {
    v.PasteMatrix(this->variables.Get_qb(), off_v, 0);
}

////

void ChBodyAddedMass::InjectVariables(ChSystemDescriptor& mdescriptor) {
    this->variables.SetDisabled(!this->IsActive());
    mdescriptor.InsertVariables(&this->variables);
}

void ChBodyAddedMass::VariablesFbReset() {
    this->variables.Get_fb().FillElem(0.0);
}

void ChBodyAddedMass::VariablesFbLoadForces(double factor) {
    // add applied forces to 'fb' vector
    this->variables.Get_fb().PasteSumVector(Xforce * factor, 0, 0);

    // add applied torques to 'fb' vector, including gyroscopic torque
    if (this->GetNoGyroTorque())
        this->variables.Get_fb().PasteSumVector((Xtorque)*factor, 3, 0);
    else
        this->variables.Get_fb().PasteSumVector((Xtorque - gyro) * factor, 3, 0);
}

void ChBodyAddedMass::VariablesFbIncrementMq() {
    this->variables.Compute_inc_Mb_v(this->variables.Get_fb(), this->variables.Get_qb());
}

void ChBodyAddedMass::VariablesQbLoadSpeed() {
    // set current speed in 'qb', it can be used by the solver when working in incremental mode
    this->variables.Get_qb().PasteVector(GetCoord_dt().pos, 0, 0);
    this->variables.Get_qb().PasteVector(GetWvel_loc(), 3, 0);
}

void ChBodyAddedMass::VariablesQbSetSpeed(double step) {
    ChCoordsys<> old_coord_dt = this->GetCoord_dt();

    // from 'qb' vector, sets body speed, and updates auxiliary data
    this->SetPos_dt(this->variables.Get_qb().ClipVector(0, 0));
    this->SetWvel_loc(this->variables.Get_qb().ClipVector(3, 0));

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

    ChVector<> newspeed = variables.Get_qb().ClipVector(0, 0);
    ChVector<> newwel = variables.Get_qb().ClipVector(3, 0);

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
    R.PasteSumVector(Xforce * c, off, 0);
    // add applied torques to 'fb' vector, including gyroscopic torque
    if (this->GetNoGyroTorque())
        R.PasteSumVector((Xtorque)*c, off + 3, 0);
    else
        R.PasteSumVector((Xtorque - gyro) * c, off + 3, 0);
}

void ChBodyAddedMass::IntLoadResidual_Mv(const unsigned int off,      // offset in R residual
                                ChVectorDynamic<>& R,        // result: the R residual, R += c*M*v
                                const ChVectorDynamic<>& w,  // the w vector
                                const double c               // a scaling factor
                                ) {
  ChMatrixDynamic<> ww = ChMatrixDynamic<>(6, 1);
  for (int i=0; i < 6; i++) {
    ww.SetElement(i, 0, w(off+i));
  }
  R.PasteSumMatrix(variables.GetMfullmass()*ww*c, off, 0);
}
}  // end namespace chrono

chrono::ChBodyAddedMass * newChBodyAddedMass()
{
  return new chrono::ChBodyAddedMass();
};
