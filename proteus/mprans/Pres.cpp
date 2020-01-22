#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "Pres.h"

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
using proteus::cppPres_base;

PYBIND11_MODULE(cPres, m)
{
    xt::import_numpy();

    py::class_<cppPres_base>(m, "Pres")
        .def(py::init(&proteus::newPres))
        .def("calculateResidual", &cppPres_base::calculateResidual)
        .def("calculateJacobian", &cppPres_base::calculateJacobian);
}
