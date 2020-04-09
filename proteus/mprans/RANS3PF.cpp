#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "RANS3PF.h"
#include "RANS3PF2D.h"

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
using proteus::cppRANS3PF2D_base;

PYBIND11_MODULE(RANS3PF, m)
{
    xt::import_numpy();

    py::class_<cppRANS3PF_base>(m, "cppRANS3PF_base")
        .def(py::init(&proteus::newRANS3PF))
        .def("calculateResidual", &cppRANS3PF_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PF_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PF_base::calculateVelocityAverage)
        .def("getBoundaryDOFs", &cppRANS3PF_base::getBoundaryDOFs);
}

PYBIND11_MODULE(RANS3PF2D, m)
{
    xt::import_numpy();

    py::class_<cppRANS3PF2D_base>(m, "cppRANS3PF2D_base")
        .def(py::init(&proteus::newRANS3PF2D))
        .def("calculateResidual", &cppRANS3PF2D_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PF2D_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PF2D_base::calculateVelocityAverage)
        .def("getBoundaryDOFs", &cppRANS3PF2D_base::getBoundaryDOFs);
}

