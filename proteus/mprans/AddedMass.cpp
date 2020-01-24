#define FORCE_IMPORT_ARRAY
#include "AddedMass.h"

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

PYBIND11_MODULE(cAddedMass, m)
{
    using proteus::cppAddedMass_base;
    using proteus::newAddedMass;

    xt::import_numpy();

    py::class_<cppAddedMass_base>(m, "cppAddedMass_base")
        .def("calculateResidual", &cppAddedMass_base::calculateResidual)
        .def("calculateJacobian", &cppAddedMass_base::calculateJacobian);

    m.def("newAddedMass", newAddedMass);
}
