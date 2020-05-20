#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "Kappa2D.h"

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
using proteus::Kappa2D_base;

PYBIND11_MODULE(cKappa2D, m)
{
    xt::import_numpy();

    py::class_<Kappa2D_base>(m, "cKappa2D_base")
        .def(py::init(&proteus::newKappa2D))
        .def("calculateResidual", &Kappa2D_base::calculateResidual)
        .def("calculateJacobian", &Kappa2D_base::calculateJacobian);
}
