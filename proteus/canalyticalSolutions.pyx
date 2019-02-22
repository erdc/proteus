# A type of -*- python -*- file
import numpy as np
cimport numpy as np
cdef extern from "analyticalSolutions.h":
    extern int cPlaneCouetteFlow_u "PlaneCouetteFlow_u"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cdiffusionSin1D "diffusionSin1D"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cdiffusionSin2D "diffusionSin2D"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cdiffusionSin3D "diffusionSin3D"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cdiffusionSin1D_r "diffusionSin1D_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cdiffusionSin2D_r "diffusionSin2D_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cdiffusionSin3D_r "diffusionSin3D_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cLinearAD_DiracIC "LinearAD_DiracIC"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cLinearAD_DiracIC_advectiveVelocity "LinearAD_DiracIC_advectiveVelocity"(int * iwork, double * rwork, int nPoints, double t, double * x, double * f)
    extern int cLinearAD_DiracIC_diffusiveVelocity "LinearAD_DiracIC_diffusiveVelocity"(int * iwork, double * rwork, int nPoints, double t, double * x, double * f)
    extern int cLinearAD_DiracIC_du "LinearAD_DiracIC_du"(int * iwork, double * rwork, int nPoints, double t, double * x, double * du)
    extern int cLinearAD_DiracIC_totalVelocity "LinearAD_DiracIC_totalVelocity"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cLinearAD_SteadyState "LinearAD_SteadyState"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cLinearADR_Decay_DiracIC "LinearADR_Decay_DiracIC"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cLinearADR_Decay_DiracIC_dr "LinearADR_Decay_DiracIC_dr"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * dr)
    extern int cLinearADR_Decay_DiracIC_r "LinearADR_Decay_DiracIC_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cLinearADR_Sine "LinearADR_Sine"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cLinearADR_Sine_advectiveVelocity "LinearADR_Sine_advectiveVelocity"(int * iwork, double * rwork, int nPoints, double t, double * x, double * f)
    extern int cLinearADR_Sine_diffusiveVelocity "LinearADR_Sine_diffusiveVelocity"(int * iwork, double * rwork, int nPoints, double t, double * x, double * f)
    extern int cLinearADR_Sine_dr "LinearADR_Sine_dr"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * dr)
    extern int cLinearADR_Sine_du "LinearADR_Sine_du"(int * iwork, double * rwork, int nPoints, double t, double * x, double * du)
    extern int cLinearADR_Sine_r "LinearADR_Sine_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cLinearADR_Sine_totalVelocity "LinearADR_Sine_totalVelocity"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cNonlinearAD_SteadyState "NonlinearAD_SteadyState"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cNonlinearADR_Decay_DiracIC "NonlinearADR_Decay_DiracIC"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cNonlinearADR_Decay_DiracIC_dr "NonlinearADR_Decay_DiracIC_dr"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * dr)
    extern int cNonlinearADR_Decay_DiracIC_r "NonlinearADR_Decay_DiracIC_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cNonlinearDAE "NonlinearDAE"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cNonlinearDAE_f "NonlinearDAE_f"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cPlanePoiseuilleFlow_u "PlanePoiseuilleFlow_u"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cPoiseuillePipeFlow "PoiseuillePipeFlow"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cPoiseuillePipeFlow_P "PoiseuillePipeFlow_P"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cpoissonsEquationExp1D "poissonsEquationExp1D"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cpoissonsEquationExp2D "poissonsEquationExp2D"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cpoissonsEquationExp3D "poissonsEquationExp3D"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cpoissonsEquationExp3D_dr "poissonsEquationExp3D_dr"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * dr)
    extern int cpoissonsEquationExp1D_r "poissonsEquationExp1D_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cpoissonsEquationExp2D_r "poissonsEquationExp2D_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cpoissonsEquationExp3D_r "poissonsEquationExp3D_r"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u, double * r)
    extern int cSTflowSphere_P "STflowSphere_P"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cSTflowSphere_Vx "STflowSphere_Vx"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cSTflowSphere_Vy "STflowSphere_Vy"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)
    extern int cSTflowSphere_Vz "STflowSphere_Vz"(int * iwork, double * rwork, int nPoints, double t, double * x, double * u)


def PlaneCouetteFlow_u(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cPlaneCouetteFlow_u( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def diffusionSin1D(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cdiffusionSin1D( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def diffusionSin2D(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cdiffusionSin2D( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def diffusionSin3D(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cdiffusionSin3D( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def diffusionSin1D_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cdiffusionSin1D_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def diffusionSin2D_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cdiffusionSin2D_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def diffusionSin3D_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cdiffusionSin3D_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def LinearAD_DiracIC(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cLinearAD_DiracIC( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def LinearAD_DiracIC_advectiveVelocity(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray f) -> int:
    return cLinearAD_DiracIC_advectiveVelocity( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > f.data)


def LinearAD_DiracIC_diffusiveVelocity(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray f) -> int:
    return cLinearAD_DiracIC_diffusiveVelocity( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > f.data)


def LinearAD_DiracIC_du(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray du) -> int:
    return cLinearAD_DiracIC_du( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > du.data)


def LinearAD_DiracIC_totalVelocity(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cLinearAD_DiracIC_totalVelocity( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def LinearAD_SteadyState(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cLinearAD_SteadyState( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def LinearADR_Decay_DiracIC(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cLinearADR_Decay_DiracIC( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def LinearADR_Decay_DiracIC_dr(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray dr) -> int:
    return cLinearADR_Decay_DiracIC_dr( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > dr.data)


def LinearADR_Decay_DiracIC_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cLinearADR_Decay_DiracIC_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def LinearADR_Sine(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cLinearADR_Sine( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def LinearADR_Sine_advectiveVelocity(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray f) -> int:
    return cLinearADR_Sine_advectiveVelocity( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > f.data)


def LinearADR_Sine_diffusiveVelocity(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray f) -> int:
    return cLinearADR_Sine_diffusiveVelocity( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > f.data)


def LinearADR_Sine_dr(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray dr) -> int:
    return cLinearADR_Sine_dr( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > dr.data)


def LinearADR_Sine_du(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray du) -> int:
    return cLinearADR_Sine_du( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > du.data)


def LinearADR_Sine_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cLinearADR_Sine_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def LinearADR_Sine_totalVelocity(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cLinearADR_Sine_totalVelocity( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def NonlinearAD_SteadyState(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cNonlinearAD_SteadyState( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def NonlinearADR_Decay_DiracIC(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cNonlinearADR_Decay_DiracIC( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def NonlinearADR_Decay_DiracIC_dr(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray dr) -> int:
    return cNonlinearADR_Decay_DiracIC_dr( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > dr.data)


def NonlinearADR_Decay_DiracIC_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cNonlinearADR_Decay_DiracIC_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def NonlinearDAE(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cNonlinearDAE( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def NonlinearDAE_f(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cNonlinearDAE_f( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def PlanePoiseuilleFlow_u(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cPlanePoiseuilleFlow_u( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def PoiseuillePipeFlow(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cPoiseuillePipeFlow( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def PoiseuillePipeFlow_P(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cPoiseuillePipeFlow_P( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def poissonsEquationExp1D(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cpoissonsEquationExp1D( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def poissonsEquationExp2D(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cpoissonsEquationExp2D( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def poissonsEquationExp3D(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cpoissonsEquationExp3D( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def poissonsEquationExp3D_dr(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray dr) -> int:
    return cpoissonsEquationExp3D_dr( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > dr.data)


def poissonsEquationExp1D_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cpoissonsEquationExp1D_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def poissonsEquationExp2D_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cpoissonsEquationExp2D_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def poissonsEquationExp3D_r(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u, np.ndarray r) -> int:
    return cpoissonsEquationExp3D_r( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data, < double * > r.data)


def STflowSphere_P(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cSTflowSphere_P( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def STflowSphere_Vx(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cSTflowSphere_Vx( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def STflowSphere_Vy(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cSTflowSphere_Vy( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)


def STflowSphere_Vz(np.ndarray iwork, np.ndarray rwork, double t, np.ndarray x, np.ndarray u) -> int:
    return cSTflowSphere_Vz( < int * > iwork.data, < double * > rwork.data, x.size/x.shape[x.ndim-1], t, < double * > x.data, < double * > u.data)
