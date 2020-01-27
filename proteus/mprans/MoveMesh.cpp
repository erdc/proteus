#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "MoveMesh.h"

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
using proteus::MoveMesh_base;

PYBIND11_MODULE(cMoveMesh, m)
{
    xt::import_numpy();

    py::class_<MoveMesh_base>(m, "cMoveMesh_base")
        .def(py::init(&proteus::newMoveMesh))
        .def("calculateResidual", &MoveMesh_base::calculateResidual)
        .def("calculateJacobian", &MoveMesh_base::calculateJacobian);
}
