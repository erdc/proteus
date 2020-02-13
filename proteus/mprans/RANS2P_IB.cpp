#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS2P_IB.h"

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
using proteus::RANS2P_IB_base;

PYBIND11_MODULE(cRANS2P_IB, m)
{
    xt::import_numpy();

    py::class_<RANS2P_IB_base>(m, "cRANS2P_IB_base")
        .def(py::init(&proteus::newRANS2P_IB))
        .def("calculateResidual"       , &RANS2P_IB_base::calculateResidual       )
        .def("calculateBeams"          , &RANS2P_IB_base::calculateBeams          )
        .def("calculateJacobian"       , &RANS2P_IB_base::calculateJacobian       )
        .def("calculateForce"          , &RANS2P_IB_base::calculateForce          )
        .def("calculateVelocityAverage", &RANS2P_IB_base::calculateVelocityAverage);
}
