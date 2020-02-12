#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS3PSed.h"
#include "RANS3PSed2D.h"

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
using proteus::cppRANS3PSed_base;
using proteus::cppRANS3PSed2D_base;

PYBIND11_MODULE(RANS3PSed, m)
{
    xt::import_numpy();

    py::class_<cppRANS3PSed_base>(m, "RANS3PSed")
        .def(py::init(&proteus::newRANS3PSed))
        .def("calculateResidual", &cppRANS3PSed_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PSed_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PSed_base::calculateVelocityAverage);
}

PYBIND11_MODULE(RANS3PSed2D, m)
{
    xt::import_numpy();

    py::class_<cppRANS3PSed2D_base>(m, "RANS3PSed2D")
        .def(py::init(&proteus::newRANS3PSed2D))
        .def("calculateResidual", &cppRANS3PSed2D_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PSed2D_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PSed2D_base::calculateVelocityAverage);
}
