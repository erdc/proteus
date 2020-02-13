#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "Kappa.h"

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
using proteus::Kappa_base;

PYBIND11_MODULE(cKappa, m)
{
    xt::import_numpy();

    py::class_<Kappa_base>(m, "cKappa_base")
        .def(py::init(&proteus::newKappa))
        .def("calculateResidual", &Kappa_base::calculateResidual)
        .def("calculateJacobian", &Kappa_base::calculateJacobian);
}
