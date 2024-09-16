// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Tristan de Lataillade, Chris Kees
// =============================================================================

#ifndef CHVARIABLESBODYADDEDMASS_H
#define CHVARIABLESBODYADDEDMASS_H

#include "chrono/solver/ChVariablesBodyOwnMass.h"

namespace chrono {

  /// Specialized class for representing a 6-DOF item for a
  /// system, that is a 3D rigid body, with mass matrix and
  /// associate variables (a 6 element vector, ex.speed)
  /// A full 6x6 matrix is used for the mass matrix, but forces
  /// still computed using the "mass" variable of ChVariablesBodyOwnMass

  class ChApi ChVariablesBodyAddedMass : public ChVariablesBodyOwnMass {

  private:
    int ndof;
    /// 6x6 mass matrix (mass and moment of inertia)
    ChMatrixDynamic<> Mmass;
    /// 6x6 added mass matrix
    ChMatrixDynamic<> Maddedmass;
    /// 6x6 full matrix (mass + added mass)
    ChMatrixDynamic<> Mfullmass;
    ChMatrixDynamic<> inv_Mfullmass;

  public:
    ChVariablesBodyAddedMass();
    virtual ~ChVariablesBodyAddedMass() {}

    /// Assignment operator: copy from other object
    ChVariablesBodyAddedMass& operator=(const ChVariablesBodyAddedMass& other);

    /* /// Access the inertia matrix */
    /* ChMatrix<>& GetMaddedmass() { return Maddedmass; } */

    /* /// Access the inertia matrix */
    /* ChMatrix<>& GetMmass() { return Mmass; } */

    /// Set the inertia matrix
    void SetMfullmass(ChMatrixDynamic<>& Mfullmass_in);

    /// Access the inertia matrix
    ChMatrixDynamic<>& GetMfullmass() { return Mfullmass; }

    /// Access the inverted inertia matrix
    ChMatrixDynamic<>& GetInvMfullmass() { return inv_Mfullmass; }

    /* /// Set the inertia matrix */
    /* void SetBodyInertia(const ChMatrix33<>& minertia); */

    /* /// Set the mass associated with translation of body */
    /* void SetBodyMass(const double mmass); */

    /* /// Get the mass associated with translation of body */
    /* virtual double GetBodyMass() const override { return ChVariablesBodyOwnMass::mmass; } */

    /* /// Access the 3x3 inertia matrix */
    /* virtual ChMatrix33<>& GetBodyInertia() override { return ChVariablesBodyOwnInertia::inertia; } */

    /* /// Set the mass associated with translation of body */
    /* void SetBodyAddedMass(ChMatrixDynamic<>& Maddedmass_in); */
    
    /// Computes the product of the inverse mass matrix by a
    /// vector, and set in result: result = [invMb]*vect
    //virtual void Compute_invMb_v(ChVectorRef result, const ChVectorConstRef vect) const override;
    virtual void ComputeMassInverseTimesVector(ChVectorRef result, ChVectorConstRef vect) const override;

    /// Computes the product of the inverse mass matrix by a
    /// vector, and increment result: result += [invMb]*vect
    virtual void Compute_inc_invMb_v(ChVectorRef result, ChVectorConstRef vect) const;

    /// Computes the product of the mass matrix by a
    /// vector, and set in result: result = [Mb]*vect
    //virtual void Compute_inc_Mb_v(ChVectorRef result, const ChVectorConstRef vect) const override;
    virtual void AddMassTimesVector(ChVectorRef results, ChVectorConstRef vect) const override;

    /// Computes the product of the corresponding block in the
    /// system matrix (ie. the mass matrix) by 'vect', scale by c_a, and add to 'result'.
    /// NOTE: the 'vect' and 'result' vectors must already have
    /// the size of the total variables&constraints in the system; the procedure
    /// will use the ChVariable offsets (that must be already updated) to know the
    /// indexes in result and vect.
    //virtual void MultiplyAndAdd(ChVectorRef result,
    //	const ChVectorConstRef vect,
    //	const double c_a) const override;
    virtual void AddMassTimesVectorInto(ChVectorRef result, ChVectorConstRef vect, const double ca) const override;

    /// Add the diagonal of the mass matrix scaled by c_a, to 'result'.
    /// NOTE: the 'result' vector must already have the size of system unknowns, ie
    /// the size of the total variables&constraints in the system; the procedure
    /// will use the ChVariable offset (that must be already updated) as index.
    //virtual void DiagonalAdd(ChVectorRef result, const double c_a) const override;
    virtual void AddMassDiagonalInto(ChVectorRef result, const double c_a) const override;

    /// Build the mass matrix (for these variables) scaled by c_a, storing
    /// it in 'storage' sparse matrix, at given column/row offset.
    /// Note, most iterative solvers don't need to know mass matrix explicitly.
    /// Optimized: doesn't fill unneeded elements except mass and 3x3 inertia.
    //virtual void Build_M(ChSparseMatrix& storage, int insrow, int inscol, const double c_a) override;
    virtual void PasteMassInto(ChSparseMatrix& mat, unsigned int insrow, unsigned int inscol, const double c_a) const override;

  };


}  // end namespace chrono

#endif
