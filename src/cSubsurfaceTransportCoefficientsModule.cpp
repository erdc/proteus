#include "Python.h"
#include "numpy/arrayobject.h"
#include "SubsurfaceTransportCoefficients.h"
#include <iostream>
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

extern "C"
{
  static PyObject *SubSurfTransCoef_calculateRusanovFluxSaturationEquationIncomp_PWC(PyObject* self, PyObject* args)
  {
    double safetyFactor,muw,mun,b;
    int failed(0);
    int nSpace,pskModelFlag,nParams,
      nElements_global,
      nElementBoundaries_element,
      nInteriorElementBoundaries_global,
      nExteriorElementBoundaries_global,
      nQuadraturePoints_element, 
      nQuadraturePoints_elementBoundary;

    PyObject *rowptr,*colind,
      *materialTypes,
      *omega,
      *Kbar,
      *rwork_psk,
      *rwork_psk_tol,
      *rwork_density_w,
      *rwork_density_n,
      *g,
      *ebq_global_qt,
      *q_lambda_bar,
      *interiorElementBoundaries,
      *exteriorElementBoundaries,
      *elementBoundaryElements,
      *elementBoundaryLocalElementBoundaries,
      *n,
      *q_u,
      *q_dV,
      *isDOFBoundary,
      *bc_u,
      *flux;
    if(!PyArg_ParseTuple(args,
			 "diiiOOOddOOdOOOOOOOiiiiiiOOOOOOOOOO",
			 &safetyFactor,
			 &nSpace,
			 &pskModelFlag,
			 &nParams,
			 &rowptr,
			 &colind,
			 &materialTypes,
			 &muw,
			 &mun,
			 &omega,
			 &Kbar,
			 &b,
			 &rwork_psk,
			 &rwork_psk_tol,
			 &rwork_density_w,
			 &rwork_density_n,
			 &g,
			 &ebq_global_qt,
			 &q_lambda_bar,
			 &nElements_global,
			 &nElementBoundaries_element,
			 &nInteriorElementBoundaries_global,
			 &nExteriorElementBoundaries_global,
			 &nQuadraturePoints_element, 
			 &nQuadraturePoints_elementBoundary,
			 &interiorElementBoundaries,
			 &exteriorElementBoundaries,
			 &elementBoundaryElements,
			 &elementBoundaryLocalElementBoundaries,
			 &n,
			 &q_u,
			 &q_dV,
			 &isDOFBoundary,
			 &bc_u,
			 &flux))
      return NULL;
    failed = calculateRusanovFluxSaturationEquationIncomp_PWC(safetyFactor,
							      nSpace,
							      pskModelFlag,
							      nParams,
							      IDATA(rowptr),
							      IDATA(colind),
							      IDATA(materialTypes),
							      muw,
							      mun,
							      DDATA(omega),
							      DDATA(Kbar),
							      b,
							      DDATA(rwork_psk),
							      DDATA(rwork_psk_tol),
							      DDATA(rwork_density_w),
							      DDATA(rwork_density_n),
							      DDATA(g),
							      DDATA(ebq_global_qt),
							      DDATA(q_lambda_bar),
							      nElements_global,
							      nElementBoundaries_element,
							      nInteriorElementBoundaries_global,
							      nExteriorElementBoundaries_global,
							      nQuadraturePoints_element, 
							      nQuadraturePoints_elementBoundary,
							      IDATA(interiorElementBoundaries),
							      IDATA(exteriorElementBoundaries),
							      IDATA(elementBoundaryElements),
							      IDATA( elementBoundaryLocalElementBoundaries),
							      DDATA(n),
							      DDATA(q_u),
							      DDATA(q_dV),
							      IDATA(isDOFBoundary),
							      DDATA(bc_u),
							      DDATA(flux));
    
    Py_INCREF(Py_None);
    return Py_None;
  }

  static PyObject* SubSurfTransCoef_piecewiseLinearTableLookup(PyObject* self, PyObject* args)

  {
    int start=0;
    double x,y,dy;
    PyObject *xv,*yv;
    if(!PyArg_ParseTuple(args,"dOO|i",
			 &x,
			 &xv,
			 &yv,
			 &start))
      return NULL;
    piecewiseLinearTableLookup(x,
			       SHAPE(xv)[0],
			       &start,
			       &y,
			       &dy,
			       DDATA(xv),
			       DDATA(yv));
    return Py_BuildValue("ddi",y,dy,start);
  }

  static PyMethodDef cSubsurfaceTransportCoefficientsMethods[] = {
    {"calculateRusanovFluxSaturationEquationIncomp_PWC",
     SubSurfTransCoef_calculateRusanovFluxSaturationEquationIncomp_PWC,
     METH_VARARGS,
     ""},
    { "piecewiseLinearTableLookup", 
      SubSurfTransCoef_piecewiseLinearTableLookup, 
      METH_VARARGS, 
      "linear table lookup routine for scalar values"},
    {NULL, NULL, 0, NULL}
  };
  PyMODINIT_FUNC initcSubsurfaceTransportCoefficients(void)
  {
    PyObject *m, *d;
    m = Py_InitModule("cSubsurfaceTransportCoefficients", cSubsurfaceTransportCoefficientsMethods);
    d = PyModule_GetDict(m);
    import_array();
  }
}//extern C
