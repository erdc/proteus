#define FORCE_IMPORT_ARRAY
#include "VOF3P.h"

namespace py = pybind11;

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

PYBIND11_MODULE(cVOF3P, m)
{
    using proteus::cppVOF3P_base;
    using proteus::newVOF3P;

    xt::import_numpy();

    py::class_<cppVOF3P_base>(m, "cppVOF3P_base")
        .def("calculateResidualElementBased", &cppVOF3P_base::calculateResidualElementBased)
        .def("calculateJacobian", &cppVOF3P_base::calculateJacobian)
        .def("FCTStep", &cppVOF3P_base::FCTStep)
        .def("calculateResidualEdgeBased", &cppVOF3P_base::calculateResidualEdgeBased);

    m.def("newVOF3P", newVOF3P);
}

