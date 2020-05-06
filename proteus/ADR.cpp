#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "ADR.h"

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
using proteus::cADR_base;

PYBIND11_MODULE(cADR, m)
{
    xt::import_numpy();

    py::class_<cppADR_base>(m, "cADR_base")
        .def(py::init(&proteus::newADR))
        .def("calculateResidual", &cppADR_base::calculateResidual)
        .def("calculateJacobian", &cppADR_base::calculateJacobian);
}
