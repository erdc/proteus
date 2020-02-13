#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS2P2D.h"

#if defined(__GNUC__) && !defined(__clang__)
    namespace workaround
    {
        inline void define_allocators()
        {
            std::allocator<int> a0;
            std::allocator<double> a1;
        }
    }
#endif

namespace py = pybind11;
using proteus::RANS2P2D_base;

PYBIND11_MODULE(cRANS2P2D, m)
{
    xt::import_numpy();

    py::class_<RANS2P2D_base>(m, "cRANS2P2D_base")
        .def(py::init(&proteus::newRANS2P2D))
        .def("calculateResidual"                    , &RANS2P2D_base::calculateResidual                     )
        .def("calculateJacobian"                    , &RANS2P2D_base::calculateJacobian                     )
        .def("calculateVelocityAverage"             , &RANS2P2D_base::calculateVelocityAverage              )
        .def("getTwoPhaseAdvectionOperator"         , &RANS2P2D_base::getTwoPhaseAdvectionOperator          )
        .def("getTwoPhaseInvScaledLaplaceOperator"  , &RANS2P2D_base::getTwoPhaseInvScaledLaplaceOperator   )
        .def("getTwoPhaseScaledMassOperator"        , &RANS2P2D_base::getTwoPhaseScaledMassOperator         );
}
