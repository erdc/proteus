#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "SW2DCV.h"

#if defined(__GNUC__) && !defined(__clang__)
namespace workaround {
inline void define_allocators() {
  std::allocator<int> a0;
  std::allocator<double> a1;
}
} // namespace workaround
#endif

namespace py = pybind11;
using proteus::SW2DCV_base;
using pybind11::return_value_policy;

PYBIND11_MODULE(cSW2DCV, m) {
  xt::import_numpy();

  py::class_<SW2DCV_base>(m, "cSW2DCV_base")
      .def(py::init(&proteus::newSW2DCV))
      .def("convexLimiting", &SW2DCV_base::convexLimiting)
      .def("calculateEdgeBasedCFL", &SW2DCV_base::calculateEdgeBasedCFL,
           return_value_policy::take_ownership)
      .def("calculateEV", &SW2DCV_base::calculateEV)
      .def("calculateResidual", &SW2DCV_base::calculateResidual)
      .def("calculateMassMatrix", &SW2DCV_base::calculateMassMatrix)
      .def("calculateLumpedMassMatrix",
           &SW2DCV_base::calculateLumpedMassMatrix);
}
