#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "MoveMesh2D.h"

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
using proteus::MoveMesh2D_base;

PYBIND11_MODULE(cMoveMesh2D, m)
{
    xt::import_numpy();

    py::class_<MoveMesh2D_base>(m, "cMoveMesh2D_base")
        .def(py::init(&proteus::newMoveMesh2D))
        .def("calculateResidual", &MoveMesh2D_base::calculateResidual)
        .def("calculateJacobian", &MoveMesh2D_base::calculateJacobian);
}
