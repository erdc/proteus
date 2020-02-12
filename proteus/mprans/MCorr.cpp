#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "MCorr.h"

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
using proteus::MCorr_base;
using pybind11::return_value_policy;

PYBIND11_MODULE(cMCorr, m)
{
    xt::import_numpy();

    py::class_<MCorr_base>(m, "cMCorr_base")
        .def(py::init(&proteus::newMCorr))
        .def("calculateResidual"                                , &MCorr_base::calculateResidual                                )
        .def("calculateJacobian"                                , &MCorr_base::calculateJacobian                                )
        .def("elementSolve"                                     , &MCorr_base::elementSolve                                     )
        .def("elementConstantSolve"                             , &MCorr_base::elementConstantSolve                             )
        .def("globalConstantRJ"                                 , &MCorr_base::globalConstantRJ, return_value_policy::take_ownership)
        .def("calculateMass"                                    , &MCorr_base::calculateMass, return_value_policy::take_ownership)
        .def("setMassQuadrature"                                , &MCorr_base::setMassQuadrature                                )
        .def("FCTStep"                                          , &MCorr_base::FCTStep                                          )
        .def("calculateMassMatrix"                              , &MCorr_base::calculateMassMatrix                              )
        .def("setMassQuadratureEdgeBasedStabilizationMethods"   , &MCorr_base::setMassQuadratureEdgeBasedStabilizationMethods   );
}
