#include "Python.h"
#include "pskRelations.h"
#include <iostream>

extern "C"
{

  static PyObject *pskRelation_vgm_calc(PyObject *self, PyObject *args)
  {
    double Sw,sw_min,sw_max,alpha,m,ns_del,eps_small;
    double rwork[4],rwork_tol[2];

    if(!PyArg_ParseTuple(args,
			 "ddddddd",
			 &Sw,
			 &sw_min,
			 &sw_max,
			 &alpha,
			 &m,
			 &eps_small,
			 &ns_del))
      return NULL;
    
    rwork[0]=sw_min; rwork[1]=sw_max; rwork[2]=alpha; rwork[3]=m;
    rwork_tol[0]=eps_small; rwork_tol[1]=ns_del;

    VGM psk(rwork);
    psk.setTolerances(rwork_tol);

    psk.calc(Sw);

    return Py_BuildValue("dddddd",
			 psk.psic,psk.krw,psk.krn,
			 psk.dpsic,psk.dkrw,psk.dkrn);

  }

  static PyObject *pskRelation_vgm_calc_from_psic(PyObject *self, PyObject *args)
  {
    double psic,sw_min,sw_max,alpha,m,ns_del,eps_small;
    double rwork[4],rwork_tol[2];

    if(!PyArg_ParseTuple(args,
			 "ddddddd",
			 &psic,
			 &sw_min,
			 &sw_max,
			 &alpha,
			 &m,
			 &eps_small,
			 &ns_del))
      return NULL;
    
    rwork[0]=sw_min; rwork[1]=sw_max; rwork[2]=alpha; rwork[3]=m;
    rwork_tol[0]=eps_small; rwork_tol[1]=ns_del;

    VGM psk(rwork);
    psk.setTolerances(rwork_tol);

    psk.calc_from_psic(psic);

    return Py_BuildValue("dddddd",
			 psk.Se,psk.krw,psk.krn,
			 psk.dSe_dpsic,psk.dkrw,psk.dkrn);

  }

  static PyMethodDef cpskRelationsMethods[] = {
    {"vgm_calc_from_sw",
     pskRelation_vgm_calc,
     METH_VARARGS,
     "pass in wetting phase saturation, Sw and sw_min,sw_max,alpha,m,eps_small,ns_del returns psic,krw,krn,dpsic,dkrw,dkrn"},
    {"vgm_calc_from_psic",
     pskRelation_vgm_calc_from_psic,
     METH_VARARGS,
     "pass in capillary pressure, Sw and sw_min,sw_max,alpha,m,eps_small,ns_del returns Se,krw,krn,dSe_psic,dkrw,dkrn"},
    {NULL,NULL, 0,NULL}
  };
  PyMODINIT_FUNC initcpskRelations(void)
  {
    PyObject *m, *d;
    m = Py_InitModule("cpskRelations", cpskRelationsMethods);
    d = PyModule_GetDict(m);
    
  }
}//end extern "C"
