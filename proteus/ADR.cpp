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

PYBIND11_MODULE(ADR, m)
{
    xt::import_numpy();

    py::class_<proteus::cppADR_base>(m, "ADR")
        .def(py::init(&proteus::newADR))
        .def("calculateResidual", &proteus::cppADR_base::calculateResidual)
        .def("calculateJacobian", &proteus::cppADR_base::calculateJacobian);
}
