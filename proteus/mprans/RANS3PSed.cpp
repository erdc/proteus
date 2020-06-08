#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS3PSed.h"

namespace py = pybind11;
using proteus::cppRANS3PSed_base;

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

PYBIND11_MODULE(cRANS3PSed, m)
{
    xt::import_numpy();

    py::class_<cppRANS3PSed_base>(m, "cppRANS3PSed_base")
        .def(py::init(&proteus::newRANS3PSed))
        .def("calculateResidual", &cppRANS3PSed_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PSed_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PSed_base::calculateVelocityAverage);
}
