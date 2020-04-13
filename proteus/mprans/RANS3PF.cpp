#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS3PF.h"

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
using proteus::cppRANS3PF_base;

PYBIND11_MODULE(cRANS3PF, m)
{
    xt::import_numpy();

    py::class_<cppRANS3PF_base>(m, "cppRANS3PF_base")
        .def(py::init(&proteus::newRANS3PF))
        .def("calculateResidual", &cppRANS3PF_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PF_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PF_base::calculateVelocityAverage)
        .def("getBoundaryDOFs", &cppRANS3PF_base::getBoundaryDOFs);

}

