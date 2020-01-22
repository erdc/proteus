#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "ElastoPlastic.h"

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
using proteus::ElastoPlastic_base;

PYBIND11_MODULE(cElastoPlastic, m)
{
    xt::import_numpy();

    py::class_<ElastoPlastic_base>(m, "cElastoPlastic_base")
        .def(py::init(&proteus::newElastoPlastic))
        .def("calculateResidual", &ElastoPlastic_base::calculateResidual)
        .def("calculateJacobian", &ElastoPlastic_base::calculateJacobian);
}
