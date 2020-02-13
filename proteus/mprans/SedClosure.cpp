#define FORCE_IMPORT_ARRAY
#include "SedClosure.h"

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

PYBIND11_MODULE(SedClosure, m)
{
    using proteus::cppHsuSedStress2D; 

    xt::import_numpy();

    py::class_<cppHsuSedStress2D>(m, "HsuSedStress")
        .def(py::init<double, double, double, double, double, double, double, double,
                      double, double, double, double, double, double, double, double,
                      double>())
        .def_readonly("aDarcy", &cppHsuSedStress2D::aDarcy_)
        .def_readonly("betaForch", &cppHsuSedStress2D::betaForch_)
        .def_readonly("grain", &cppHsuSedStress2D::grain_)
        .def_readonly("packFraction", &cppHsuSedStress2D::packFraction_)
        .def_readonly("packMargin", &cppHsuSedStress2D::packMargin_)
        .def_readonly("maxFraction", &cppHsuSedStress2D::maxFraction_)
        .def_readonly("frFraction", &cppHsuSedStress2D::frFraction_)
        .def_readonly("sigmaC", &cppHsuSedStress2D::sigmaC_)
        .def_readonly("C3e", &cppHsuSedStress2D::C3e_)
        .def_readonly("C4e", &cppHsuSedStress2D::C4e_)
        .def_readonly("eR", &cppHsuSedStress2D::eR_)
        .def_readonly("fContact", &cppHsuSedStress2D::fContact_)
        .def_readonly("mContact", &cppHsuSedStress2D::mContact_)
        .def_readonly("nContact", &cppHsuSedStress2D::nContact_)
        .def_readonly("angFriction", &cppHsuSedStress2D::angFriction_)
        .def_readonly("vos_limiter", &cppHsuSedStress2D::vos_limiter_)
        .def_readonly("mu_fr_limiter", &cppHsuSedStress2D::mu_fr_limiter_)
        .def("betaCoeff", &cppHsuSedStress2D::xt_betaCoeff)
        .def("gs0", &cppHsuSedStress2D::gs0)
        .def("kappa_sed1", &cppHsuSedStress2D::xt_kappa_sed1)
        .def("dkappa_sed1_dk", &cppHsuSedStress2D::xt_dkappa_sed1_dk)
        .def("deps_sed_deps", &cppHsuSedStress2D::xt_deps_sed_deps)
        .def("psc", &cppHsuSedStress2D::psc)
        .def("psc_term", &cppHsuSedStress2D::psc_term)
        .def("dpsc_term_dtheta", &cppHsuSedStress2D::dpsc_term_dtheta)
        .def("mu_sc", &cppHsuSedStress2D::mu_sc)
        .def("mu_fr", &cppHsuSedStress2D::mu_fr)
        .def("l_sc", &cppHsuSedStress2D::l_sc)
        .def("tausc_term_theta", &cppHsuSedStress2D::tausc_term_theta)
        .def("gamma_s", &cppHsuSedStress2D::gamma_s)
        .def("dgamma_s_dtheta", &cppHsuSedStress2D::dgamma_s_dtheta)
        .def("jint1", &cppHsuSedStress2D::xt_jint1)
        .def("jint2", &cppHsuSedStress2D::xt_jint2)
        .def("djint2_dtheta", &cppHsuSedStress2D::xt_djint2_dtheta)
        .def("k_diff", &cppHsuSedStress2D::k_diff)
        .def("p_friction", &cppHsuSedStress2D::p_friction)
        .def("gradp_friction", &cppHsuSedStress2D::gradp_friction)
        .def("mIntFluid", &cppHsuSedStress2D::xt_mIntFluid)
        .def("mIntSolid", &cppHsuSedStress2D::xt_mIntSolid)
        .def("mIntgradC", &cppHsuSedStress2D::xt_mIntgradC)
        .def("dmInt_duFluid", &cppHsuSedStress2D::xt_dmInt_duFluid)
        .def("dmInt_duSolid", &cppHsuSedStress2D::xt_dmInt_duSolid)
        .def("p_s", &cppHsuSedStress2D::p_s);
}
