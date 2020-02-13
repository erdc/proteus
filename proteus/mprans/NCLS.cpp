#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "NCLS.h"

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
using proteus::NCLS_base;

PYBIND11_MODULE(cNCLS, m)
{
    xt::import_numpy();

    py::class_<NCLS_base>(m, "cNCLS_base")
        .def(py::init(&proteus::newNCLS))
        .def("calculateResidual"                    , &NCLS_base::calculateResidual                     )
        .def("calculateJacobian"                    , &NCLS_base::calculateJacobian                     )
        .def("calculateWaterline"                   , &NCLS_base::calculateWaterline                    )
        .def("calculateRedistancingResidual"        , &NCLS_base::calculateRedistancingResidual         )
        .def("calculateRhsSmoothing"                , &NCLS_base::calculateRhsSmoothing                 )
        .def("calculateResidual_entropy_viscosity"  , &NCLS_base::calculateResidual_entropy_viscosity   )
        .def("calculateMassMatrix"                  , &NCLS_base::calculateMassMatrix                   )
        .def("calculateSmoothingMatrix"             , &NCLS_base::calculateSmoothingMatrix              );
}
