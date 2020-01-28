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

PYBIND11_MODULE(cRANS3PF, m)
{
    using proteus::cppRANS3PF_base;
    using proteus::newRANS3PF;
    using proteus::cppRANS3PF2D_base;
    using proteus::newRANS3PF2D;

    xt::import_numpy();

    py::class_<cppRANS3PF_base>(m, "RANS3PF")
        .def("calculateResidual", &cppRANS3PF_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PF_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PF_base::calculateVelocityAverage)
        .def("getBoundaryDOFs", &cppRANS3PF_base::getBoundaryDOFs);

    m.def("newRANS3PF", newRANS3PF);

    py::class_<cppRANS3PF2D_base>(m, "RANS3PF2D")
        .def("calculateResidual", &cppRANS3PF2D_base::calculateResidual)
        .def("calculateJacobian", &cppRANS3PF2D_base::calculateJacobian)
        .def("calculateVelocityAverage", &cppRANS3PF2D_base::calculateVelocityAverage)
        .def("getBoundaryDOFs", &cppRANS3PF2D_base::getBoundaryDOFs);

    m.def("newRANS3PF2D", newRANS3PF2D);
}

