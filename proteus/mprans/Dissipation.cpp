#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "Dissipation.h"

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
using proteus::Dissipation_base;

PYBIND11_MODULE(cDissipation, m)
{
    xt::import_numpy();

    py::class_<Dissipation_base>(m, "cDissipation_base")
        .def(py::init(&proteus::newDissipation))
        .def("calculateResidual", &Dissipation_base::calculateResidual)
        .def("calculateJacobian", &Dissipation_base::calculateJacobian);
}
