#include <cmath>
#include <iostream>
#include "CompKernel.h"

namespace proteus
{
  class RANS2P_Op_Builder
  {
  public:
    void attachTwoPhaseInvScaledMassOperator ();
  };

  void
  RANS2P_Op_Builder::attachTwoPhaseInvScaledMassOperator()
  {
    std::cout << "Function called\n";
  }
}
