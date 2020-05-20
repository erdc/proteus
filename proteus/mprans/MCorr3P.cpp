#define FORCE_IMPORT_ARRAY
#include "MCorr3P.h"

namespace xt
{
#if defined(__GNUC__) && !defined(__clang__)
    namespace workaround
    {
        inline void complex_allocator()
        {
           std::allocator<int> ai;
           std::allocator<double> ad;
        }
    }
#endif
}

namespace py = pybind11;

PYBIND11_MODULE(cMCorr3P, m)
{
    using proteus::cppMCorr3P_base;
    using proteus::newMCorr3P;

    xt::import_numpy();

    py::class_<cppMCorr3P_base>(m, "cppMCorr3P_base")
        .def("calculateResidual", &cppMCorr3P_base::calculateResidual)
        .def("calculateJacobian", &cppMCorr3P_base::calculateJacobian)
        .def("elementSolve", &cppMCorr3P_base::elementSolve)
        .def("elementConstantSolve", &cppMCorr3P_base::elementConstantSolve)
        .def("globalConstantRJ", &cppMCorr3P_base::globalConstantRJ)
        .def("calculateMass", &cppMCorr3P_base::calculateMass)
        .def("setMassQuadrature", &cppMCorr3P_base::setMassQuadrature)
        .def("calculateStiffnessMatrix", &cppMCorr3P_base::calculateStiffnessMatrix);

    m.def("newMCorr3P", newMCorr3P);
}
