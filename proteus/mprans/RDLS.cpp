#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RDLS.h"

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
using proteus::RDLS_base;
using pybind11::return_value_policy;

PYBIND11_MODULE(cRDLS, m)
{
    xt::import_numpy();

    py::class_<RDLS_base>(m, "cRDLS_base")
        .def(py::init(&proteus::newRDLS))
        .def("calculateResidual"                , &RDLS_base::calculateResidual                 )
        .def("calculateJacobian"                , &RDLS_base::calculateJacobian                 )
        .def("calculateResidual_ellipticRedist" , &RDLS_base::calculateResidual_ellipticRedist  )
        .def("calculateJacobian_ellipticRedist" , &RDLS_base::calculateJacobian_ellipticRedist  )
        .def("normalReconstruction"             , &RDLS_base::normalReconstruction              )
        .def("calculateMetricsAtEOS"            , &RDLS_base::calculateMetricsAtEOS, return_value_policy::take_ownership);
}
