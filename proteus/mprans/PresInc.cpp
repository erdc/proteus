#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "PresInc.h"

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
using proteus::cppPresInc_base;

PYBIND11_MODULE(cPresInc, m)
{
    xt::import_numpy();

    py::class_<cppPresInc_base>(m, "PresInc")
        .def(py::init(&proteus::newPresInc))
        .def("calculateResidual", &cppPresInc_base::calculateResidual)
        .def("calculateJacobian", &cppPresInc_base::calculateJacobian);
}
