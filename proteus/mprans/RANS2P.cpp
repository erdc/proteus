#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS2P.h"

namespace py = pybind11;
using proteus::RANS2P_base;

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

PYBIND11_MODULE(cRANS2P, m)
{
    xt::import_numpy();

    py::class_<RANS2P_base>(m, "cRANS2P_base")
        .def(py::init(&proteus::newRANS2P))
        .def("calculateResidual"                    , &RANS2P_base::calculateResidual                     )
        .def("calculateJacobian"                    , &RANS2P_base::calculateJacobian                     )
        .def("calculateVelocityAverage"             , &RANS2P_base::calculateVelocityAverage              )
        .def("getTwoPhaseAdvectionOperator"         , &RANS2P_base::getTwoPhaseAdvectionOperator          )
        .def("getTwoPhaseInvScaledLaplaceOperator"  , &RANS2P_base::getTwoPhaseInvScaledLaplaceOperator   )
        .def("getTwoPhaseScaledMassOperator"        , &RANS2P_base::getTwoPhaseScaledMassOperator         )
        .def("step6DOF"        , &RANS2P_base::step6DOF);
}
