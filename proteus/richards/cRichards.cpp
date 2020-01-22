#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "Richards.h"

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
using proteus::Richards_base;

PYBIND11_MODULE(cRichards, m)
{
    xt::import_numpy();

    py::class_<Richards_base>(m, "cRichards_base")
        .def(py::init(&proteus::newRichards))
        .def("calculateResidual", &Richards_base::calculateResidual)
        .def("calculateJacobian", &Richards_base::calculateJacobian);
}
