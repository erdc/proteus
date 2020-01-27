#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "VOF.h"

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
using proteus::VOF_base;

PYBIND11_MODULE(cVOF, m)
{
    xt::import_numpy();

    py::class_<VOF_base>(m, "cVOF_base")
        .def(py::init(&proteus::newVOF))
        .def("calculateResidualElementBased"    , &VOF_base::calculateResidualElementBased  )
        .def("calculateJacobian"                , &VOF_base::calculateJacobian              )
        .def("FCTStep"                          , &VOF_base::FCTStep                        )
        .def("calculateResidualEdgeBased"       , &VOF_base::calculateResidualEdgeBased     );
}
