#define FORCE_IMPORT_ARRAY
#include "NCLS3P.h"

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

PYBIND11_MODULE(cNCLS3P, m)
{
    using proteus::cppNCLS3P_base;
    using proteus::newNCLS3P;

    xt::import_numpy();

    py::class_<cppNCLS3P_base>(m, "cppNCLS3P_base")
        .def("calculateResidual", &cppNCLS3P_base::calculateResidual)
        .def("calculateJacobian", &cppNCLS3P_base::calculateJacobian)
        .def("calculateWaterline", &cppNCLS3P_base::calculateWaterline);

    m.def("newNCLS3P", newNCLS3P);
}

