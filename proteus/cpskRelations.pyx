# A type of -*- python -*- file
import numpy as np
cimport numpy as np
cimport pskRelations

def vgm_cal_from_sw(double Sw,
                    double sw_min,
                    double sw_max,
                    double alpha,
                    double m,
                    double eps_small,
                    double ns_del):
    cdef double rwork[4], rwork_tol[2]
    rwork[0]=sw_min
    rwork[1]=sw_max
    rwork[2]=alpha
    rwork[3]=m;
    rwork_tol[0]=eps_small
    rwork_tol[1]=ns_del
    cdef pskRelations.VGM psk = pskRelations.VGM(rwork)
    psk.setTolerances(rwork_tol)
    psk.calc(Sw);
    return (psk.psic,psk.krw,psk.krn,psk.dpsic,psk.dkrw,psk.dkrn)

def vgm_calc_from_psic(double psic,
		       double sw_min,
		       double sw_max,
		       double alpha,
		       double m,
		       double eps_small,
		       double ns_del):
    cdef double rwork[4],rwork_tol[2];
    rwork[0]=sw_min
    rwork[1]=sw_max
    rwork[2]=alpha
    rwork[3]=m;
    rwork_tol[0]=eps_small
    rwork_tol[1]=ns_del
    cdef pskRelations.VGM psk = pskRelations.VGM(rwork)
    psk.setTolerances(rwork_tol)
    psk.calc_from_psic(psic)
    return (psk.Se,psk.krw,psk.krn,psk.dSe_dpsic,psk.dkrw,psk.dkrn)
