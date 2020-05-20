#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "PresInit.h"

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
using proteus::cppPresInit_base;

PYBIND11_MODULE(cPresInit, m)
{
    xt::import_numpy();

    py::class_<cppPresInit_base>(m, "PresInit")
        .def(py::init(&proteus::newPresInit))
        .def("calculateResidual", &cppPresInit_base::calculateResidual)
        .def("calculateJacobian", &cppPresInit_base::calculateJacobian);
}
