#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "CLSVOF.h"

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
using proteus::CLSVOF_base;

PYBIND11_MODULE(cCLSVOF, m)
{
    xt::import_numpy();

    py::class_<CLSVOF_base>(m, "cCLSVOF_base")
        .def(py::init(&proteus::newCLSVOF))
        .def("calculateResidual"        , &CLSVOF_base::calculateResidual        )
        .def("calculateJacobian"        , &CLSVOF_base::calculateJacobian        )
        .def("calculateMetricsAtEOS"    , &CLSVOF_base::calculateMetricsAtEOS    )    
        .def("calculateMetricsAtETS"    , &CLSVOF_base::calculateMetricsAtETS    )    
        .def("normalReconstruction"     , &CLSVOF_base::normalReconstruction     )   
        .def("calculateRhsL2Proj"       , &CLSVOF_base::calculateRhsL2Proj       ) 
        .def("calculateLumpedMassMatrix", &CLSVOF_base::calculateLumpedMassMatrix)        
        .def("assembleSpinUpSystem"     , &CLSVOF_base::assembleSpinUpSystem     )   
        .def("FCTStep"                  , &CLSVOF_base::FCTStep                  );
}
