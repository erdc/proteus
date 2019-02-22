#ifndef CHBODYADDEDMASS_H
#define CHBODYADDEDMASS_H

#include "chrono/physics/ChBody.h"
#include "ChVariablesBodyAddedMass.h"

namespace chrono {
	class ChBodyAddedMass : public ChBody {
      protected:
          ChVariablesBodyAddedMass variables;
      public:
          ChBodyAddedMass();
          virtual ~ChBodyAddedMass() {}
    void SetMass(double newmass);
    void SetInertia(const ChMatrix33<>& iner);
    void SetInertiaXX(const ChVector<>& iner);
    void SetInertiaXY(const ChVector<>& iner);
    /* ChVector<> GetInertiaXX(); */
    /* ChVector<> GetInertiaXY(); */
    /* double GetMass() { return variables.GetBodyMass(); }; */
    /* const ChMatrix33<>& GetInertia() { return variables.GetBodyInertia(); } */
    void SetMfullmass(ChMatrixDynamic<> Mfullmass_in);
    void SetInvMfullmass(ChMatrixDynamic<> inv_Mfullmass_in);
    ChVariables& Variables() override { return variables; } 
    ChVariablesBodyOwnMass& VariablesBody() override { return variables; }
    ChVariablesBodyAddedMass& VariablesBodyAddedMass() { return variables; }
    //
    // STATE FUNCTIONS
    //

    // (override/implement interfaces for global state vectors, see ChPhysicsItem for comments.)

    virtual void IntToDescriptor(const unsigned int off_v,
                                 const ChStateDelta& v,
                                 const ChVectorDynamic<>& R,
                                 const unsigned int off_L,
                                 const ChVectorDynamic<>& L,
                                 const ChVectorDynamic<>& Qc) override;
    virtual void IntFromDescriptor(const unsigned int off_v,
                                   ChStateDelta& v,
                                   const unsigned int off_L,
                                   ChVectorDynamic<>& L) override;
    virtual void IntLoadResidual_F(const unsigned int off, ChVectorDynamic<>& R, const double c) override;

    virtual void IntLoadResidual_Mv(const unsigned int off,
                                    ChVectorDynamic<>& R,
                                    const ChVectorDynamic<>& w,
                                    const double c) override;

    //
    // SOLVER FUNCTIONS
    //

    // Override/implement solver system functions of ChPhysicsItem
    // (to assemble/manage data for system solver)

    /// Sets the 'fb' part of the encapsulated ChVariablesBodyOwnMass to zero.
    virtual void VariablesFbReset() override;

    /// Adds the current forces applied to body (including gyroscopic torque) in
    /// encapsulated ChVariablesBody, in the 'fb' part: qf+=forces*factor
    virtual void VariablesFbLoadForces(double factor = 1) override;

    /// Initialize the 'qb' part of the ChVariablesBody with the
    /// current value of body speeds. Note: since 'qb' is the unknown, this
    /// function seems unnecessary, unless used before VariablesFbIncrementMq()
    virtual void VariablesQbLoadSpeed() override;

    /// Adds M*q (masses multiplied current 'qb') to Fb, ex. if qb is initialized
    /// with v_old using VariablesQbLoadSpeed, this method can be used in
    /// timestepping schemes that do: M*v_new = M*v_old + forces*dt
    virtual void VariablesFbIncrementMq() override;

    /// Fetches the body speed (both linear and angular) from the
    /// 'qb' part of the ChVariablesBody (does not updates the full body&markers state)
    /// and sets it as the current body speed.
    /// If 'step' is not 0, also computes the approximate acceleration of
    /// the body using backward differences, that is  accel=(new_speed-old_speed)/step.
    /// Mostly used after the solver provided the solution in ChVariablesBody .
    virtual void VariablesQbSetSpeed(double step = 0) override;

    /// Increment body position by the 'qb' part of the ChVariablesBody,
    /// multiplied by a 'step' factor.
    ///     pos+=qb*step
    /// If qb is a speed, this behaves like a single step of 1-st order
    /// numerical integration (Eulero integration).
    /// Does not automatically update markers & forces.
    virtual void VariablesQbIncrementPosition(double step) override;

    /// Tell to a system descriptor that there are variables of type
    /// ChVariables in this object (for further passing it to a solver)
    virtual void InjectVariables(ChSystemDescriptor& mdescriptor) override;

    };

}  // end namespace chrono

chrono::ChBodyAddedMass * newChBodyAddedMass();

#endif

