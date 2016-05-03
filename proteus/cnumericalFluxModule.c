#include "Python.h"
#include "numpy/arrayobject.h"
#include "numericalFlux.h"
/** \file cnumericalFluxModule.c
    \defgroup cnumericalFlux cnumericalFlux
    \brief Python interface to numericalFlux library
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject*
cnumericalFluxCalculateInteriorNumericalAdvectiveFluxConvexOneSonicPoint(PyObject* self,
									 PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*dadvectiveFlux_right;
  double sonicPoint,sonicFlux;
  if(!PyArg_ParseTuple(args,"ddOOOOOOOOOO",
		       &sonicPoint,
		       &sonicFlux,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left,
                       &dadvectiveFlux_right))
    return NULL;
  calculateInteriorNumericalAdvectiveFluxConvexOneSonicPoint(sonicPoint,sonicFlux,
							     SHAPE(interiorElementBoundaries)[0],
							     SHAPE(f)[1],
							     SHAPE(f)[2],
							     SHAPE(f)[3],
							     IDATA(interiorElementBoundaries),
							     IDATA(elementBoundaryElements),
							     IDATA(elementBoundaryLocalElementBoundaries),
							     DDATA(n),
							     DDATA(u),
							     DDATA(f),
							     DDATA(df),
							     DDATA(advectiveFlux),
							     DDATA(dadvectiveFlux_left),
							     DDATA(dadvectiveFlux_right));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalAdvectiveFluxRusanov(PyObject* self,
							     PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*u,*f,*df,*df_element,*advectiveFlux,*dadvectiveFlux_left,*dadvectiveFlux_right;
  double safetyFactor;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOO",
		       &safetyFactor,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &f,
                       &df,
                       &df_element,
                       &advectiveFlux,
                       &dadvectiveFlux_left,
                       &dadvectiveFlux_right))
    return NULL;
  calculateInteriorNumericalAdvectiveFluxRusanov(safetyFactor,
						 SHAPE(interiorElementBoundaries)[0],
						 SHAPE(f)[1],
						 SHAPE(f)[2],
						 SHAPE(df_element)[1],
						 SHAPE(f)[3],
						 IDATA(interiorElementBoundaries),
						 IDATA(elementBoundaryElements),
						 IDATA(elementBoundaryLocalElementBoundaries),
						 DDATA(n),
						 DDATA(u),
						 DDATA(f),
						 DDATA(df),
						 DDATA(df_element),
						 DDATA(advectiveFlux),
						 DDATA(dadvectiveFlux_left),
						 DDATA(dadvectiveFlux_right));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxRusanov(PyObject* self, 
							     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*bc_u,*bc_f,*bc_df,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*inflowFlag,*isDOFBoundary,*df_element;
  double safetyFactor;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOOOOO",
		       &safetyFactor,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &inflowFlag,
                       &n,
                       &bc_u,
                       &bc_f,
                       &bc_df,
                       &u,
                       &f,
                       &df,
		       &df_element,
                       &advectiveFlux,
                       &dadvectiveFlux_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFluxRusanov(safetyFactor,
						     SHAPE(exteriorElementBoundaries)[0],
						     SHAPE(f)[1],
						     SHAPE(f)[2],
						     SHAPE(df_element)[1],
						     SHAPE(f)[3],
						     IDATA(exteriorElementBoundaries),
						     IDATA(elementBoundaryElements),
						     IDATA(elementBoundaryLocalElementBoundaries),
						     IDATA(isDOFBoundary),
						     IDATA(inflowFlag),
						     DDATA(n),
						     DDATA(bc_u),
						     DDATA(bc_f),
						     DDATA(bc_df),
						     DDATA(u),
						     DDATA(f),
						     DDATA(df),
						     DDATA(df_element),
						     DDATA(advectiveFlux),
						     DDATA(dadvectiveFlux_left));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFluxRusanov(safetyFactor,
							   SHAPE(exteriorElementBoundaries)[0],
							   SHAPE(f)[1],
							   SHAPE(df_element)[1],
							   SHAPE(f)[2],
							   IDATA(exteriorElementBoundaries),
							   IDATA(elementBoundaryElements),
							   IDATA(elementBoundaryLocalElementBoundaries),
							   IDATA(isDOFBoundary),
							   IDATA(inflowFlag),
							   DDATA(n),
							   DDATA(bc_u),
							   DDATA(bc_f),
							   DDATA(bc_df),
							   DDATA(u),
							   DDATA(f),
							   DDATA(df),
							   DDATA(df_element),
							   DDATA(advectiveFlux),
							   DDATA(dadvectiveFlux_left));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateInteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(PyObject* self,
										PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*u,*f,*lambda_bar_element,*advectiveFlux;
  double safetyFactor;
  if(!PyArg_ParseTuple(args,"dOOOOOOOO",
		       &safetyFactor,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &f,
                       &lambda_bar_element,
                       &advectiveFlux))
    return NULL;
  calculateInteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(safetyFactor,
								    SHAPE(interiorElementBoundaries)[0],
								    SHAPE(f)[1],
								    SHAPE(f)[2],
								    SHAPE(lambda_bar_element)[1],
								    SHAPE(f)[3],
								    IDATA(interiorElementBoundaries),
								    IDATA(elementBoundaryElements),
								    IDATA(elementBoundaryLocalElementBoundaries),
								    DDATA(n),
								    DDATA(u),
								    DDATA(f),
								    DDATA(lambda_bar_element),
								    DDATA(advectiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(PyObject* self, 
										PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*bc_u,*bc_f,*u,*f,*advectiveFlux,*inflowFlag,*isDOFBoundary,*lambda_bar_element;
  double safetyFactor;
  if(!PyArg_ParseTuple(args,"dOOOOOOOOOOOO",
		       &safetyFactor,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &inflowFlag,
                       &n,
                       &bc_u,
                       &bc_f,
                       &u,
                       &f,
		       &lambda_bar_element,
                       &advectiveFlux))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(safetyFactor,
						     SHAPE(exteriorElementBoundaries)[0],
						     SHAPE(f)[1],
						     SHAPE(f)[2],
						     SHAPE(lambda_bar_element)[1],
						     SHAPE(f)[3],
						     IDATA(exteriorElementBoundaries),
						     IDATA(elementBoundaryElements),
						     IDATA(elementBoundaryLocalElementBoundaries),
						     IDATA(isDOFBoundary),
						     IDATA(inflowFlag),
						     DDATA(n),
						     DDATA(bc_u),
						     DDATA(bc_f),
						     DDATA(u),
						     DDATA(f),
						     DDATA(lambda_bar_element),
									DDATA(advectiveFlux));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound(safetyFactor,
							   SHAPE(exteriorElementBoundaries)[0],
							   SHAPE(f)[1],
							   SHAPE(lambda_bar_element)[1],
							   SHAPE(f)[2],
							   IDATA(exteriorElementBoundaries),
							   IDATA(elementBoundaryElements),
							   IDATA(elementBoundaryLocalElementBoundaries),
							   IDATA(isDOFBoundary),
							   IDATA(inflowFlag),
							   DDATA(n),
							   DDATA(bc_u),
							   DDATA(bc_f),
							   DDATA(u),
							   DDATA(f),
							   DDATA(lambda_bar_element),
							   DDATA(advectiveFlux));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateInteriorLesaintRaviartNumericalFlux(PyObject* self,
							   PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*u,*H,*dH,*HamiltonJacobiFlux,*dHamiltonJacobiFlux_left,*dHamiltonJacobiFlux_right;
  int speedEvalFlag=0;
  if(!PyArg_ParseTuple(args,"iOOOOOOOOOO",
		       &speedEvalFlag,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &H,
                       &dH,
                       &HamiltonJacobiFlux,
                       &dHamiltonJacobiFlux_left,
                       &dHamiltonJacobiFlux_right))
    return NULL;
  calculateInteriorLesaintRaviartNumericalFlux(SHAPE(interiorElementBoundaries)[0],
					       SHAPE(dH)[1],
					       SHAPE(dH)[2],
					       SHAPE(dH)[3],
					       speedEvalFlag,
					       IDATA(interiorElementBoundaries),
					       IDATA(elementBoundaryElements),
					       IDATA(elementBoundaryLocalElementBoundaries),
					       DDATA(n),
					       DDATA(u),
					       DDATA(H),
					       DDATA(dH),
					       DDATA(HamiltonJacobiFlux),
					       DDATA(dHamiltonJacobiFlux_left),
					       DDATA(dHamiltonJacobiFlux_right));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateExteriorLesaintRaviartNumericalFlux(PyObject* self, 
							   PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*bc_u,*bc_H,*bc_dH,*u,*H,*dH,*HamiltonJacobiFlux,*dHamiltonJacobiFlux_left,*inflowFlag,*isDOFBoundary;
  int speedEvalFlag = 0;
  if(!PyArg_ParseTuple(args,"iOOOOOOOOOOOOOO",
		       &speedEvalFlag,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &inflowFlag,
                       &n,
                       &bc_u,
                       &bc_H,
                       &bc_dH,
                       &u,
                       &H,
                       &dH,
                       &HamiltonJacobiFlux,
                       &dHamiltonJacobiFlux_left))
    return NULL;
  if (ND(dH) > 3)
    {
      calculateExteriorLesaintRaviartNumericalFlux(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(dH)[1],
						   SHAPE(dH)[2],
						   SHAPE(dH)[3],
						   speedEvalFlag,
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary),
						   IDATA(inflowFlag),
						   DDATA(n),
						   DDATA(bc_u),
						   DDATA(bc_H),
						   DDATA(bc_dH),
						   DDATA(u),
						   DDATA(H),
						   DDATA(dH),
						   DDATA(HamiltonJacobiFlux),
						   DDATA(dHamiltonJacobiFlux_left));
    }
  else
    {
      calculateGlobalExteriorLesaintRaviartNumericalFlux(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(dH)[1],
							 SHAPE(dH)[2],
							 speedEvalFlag,
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(isDOFBoundary),
							 IDATA(inflowFlag),
							 DDATA(n),
							 DDATA(bc_u),
							 DDATA(bc_H),
							 DDATA(bc_dH),
							 DDATA(u),
							 DDATA(H),
							 DDATA(dH),
							 DDATA(HamiltonJacobiFlux),
							 DDATA(dHamiltonJacobiFlux_left));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFlux_NoBC(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*inflowFlag;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &inflowFlag,
                       &n,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFlux_NoBC(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(f)[1],
						   SHAPE(f)[2],
						   SHAPE(f)[3],
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(inflowFlag),
						   DDATA(n),
						   DDATA(f),
						   DDATA(df),
						   DDATA(advectiveFlux),
						   DDATA(dadvectiveFlux_left));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFlux_NoBC(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(f)[1],
							 SHAPE(f)[2],
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(inflowFlag),
							 DDATA(n),
							 DDATA(f),
							 DDATA(df),
							 DDATA(advectiveFlux),
							 DDATA(dadvectiveFlux_left));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFlux(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*bc_u,*bc_f,*bc_df,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*inflowFlag,*isDOFBoundary;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &inflowFlag,
                       &n,
                       &bc_u,
                       &bc_f,
                       &bc_df,
                       &u,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFlux(SHAPE(exteriorElementBoundaries)[0],
					      SHAPE(f)[1],
					      SHAPE(f)[2],
					      SHAPE(f)[3],
					      IDATA(exteriorElementBoundaries),
					      IDATA(elementBoundaryElements),
					      IDATA(elementBoundaryLocalElementBoundaries),
					      IDATA(isDOFBoundary),
					      IDATA(inflowFlag),
					      DDATA(n),
					      DDATA(bc_u),
					      DDATA(bc_f),
					      DDATA(bc_df),
					      DDATA(u),
					      DDATA(f),
					      DDATA(df),
					      DDATA(advectiveFlux),
					      DDATA(dadvectiveFlux_left));
    }
  else
    {

      calculateGlobalExteriorNumericalAdvectiveFlux(SHAPE(exteriorElementBoundaries)[0],
						    SHAPE(f)[1],
						    SHAPE(f)[2],
						    IDATA(exteriorElementBoundaries),
						    IDATA(elementBoundaryElements),
						    IDATA(elementBoundaryLocalElementBoundaries),
						    IDATA(isDOFBoundary),
						    IDATA(inflowFlag),
						    DDATA(n),
						    DDATA(bc_u),
						    DDATA(bc_f),
						    DDATA(bc_df),
						    DDATA(u),
						    DDATA(f),
						    DDATA(df),
						    DDATA(advectiveFlux),
						    DDATA(dadvectiveFlux_left));
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFlux_free(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*bc_u,*bc_f,*bc_df,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*inflowFlag,*isDOFBoundary;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &inflowFlag,
                       &n,
                       &bc_u,
                       &bc_f,
                       &bc_df,
                       &u,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFlux_free(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(f)[1],
						   SHAPE(f)[2],
						   SHAPE(f)[3],
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary),
						   IDATA(inflowFlag),
						   DDATA(n),
						   DDATA(bc_u),
						   DDATA(bc_f),
						   DDATA(bc_df),
						   DDATA(u),
						   DDATA(f),
						   DDATA(df),
						   DDATA(advectiveFlux),
						   DDATA(dadvectiveFlux_left));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFlux_free(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(f)[1],
							 SHAPE(f)[2],
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(isDOFBoundary),
							 IDATA(inflowFlag),
							 DDATA(n),
							 DDATA(bc_u),
							 DDATA(bc_f),
							 DDATA(bc_df),
							 DDATA(u),
							 DDATA(f),
							 DDATA(df),
							 DDATA(advectiveFlux),
							 DDATA(dadvectiveFlux_left));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokesP2D(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *bc_f,
    *bc_fpu,
    *bc_fpv,
    *f,
    *fpu,
    *fpv,
    *df_du,
    *df_dv,
    *dfpu_dp,
    *dfpv_dp,
    *advectiveFlux,
    *advectiveFluxpu,
    *advectiveFluxpv,
    *dadvectiveFlux_du_left,
    *dadvectiveFlux_dv_left,
    *dadvectiveFluxpu_dp_left,
    *dadvectiveFluxpv_dp_left,
    *isDOFBoundary_p,
    *isDOFBoundary_u,
    *isDOFBoundary_v;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_p,
                       &isDOFBoundary_u,
                       &isDOFBoundary_v,
                       &n,
                       &bc_f,
                       &bc_fpu,
                       &bc_fpv,
                       &f,
                       &fpu,
                       &fpv,
                       &df_du,
                       &df_dv,
                       &dfpu_dp,
                       &dfpv_dp,
                       &advectiveFlux,
                       &advectiveFluxpu,
                       &advectiveFluxpv,
                       &dadvectiveFlux_du_left,
                       &dadvectiveFlux_dv_left,
                       &dadvectiveFluxpu_dp_left,
                       &dadvectiveFluxpv_dp_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFluxStokesP2D(SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(f)[1],
						      SHAPE(f)[2],
						      SHAPE(f)[3],
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(isDOFBoundary_p),
						      IDATA(isDOFBoundary_u),
						      IDATA(isDOFBoundary_v),
						      DDATA(n),
						      DDATA(bc_f),
						      DDATA(bc_fpu),
						      DDATA(bc_fpv),
						      DDATA(f),
						      DDATA(fpu),
						      DDATA(fpv),
						      DDATA(df_du),
						      DDATA(df_dv),
						      DDATA(dfpu_dp),
						      DDATA(dfpv_dp),
						      DDATA(advectiveFlux),
						      DDATA(advectiveFluxpu),
						      DDATA(advectiveFluxpv),
						      DDATA(dadvectiveFlux_du_left),
						      DDATA(dadvectiveFlux_dv_left),
						      DDATA(dadvectiveFluxpu_dp_left),
						      DDATA(dadvectiveFluxpv_dp_left));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFluxStokesP2D(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(f)[1],
							    SHAPE(f)[2],
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary_p),
							    IDATA(isDOFBoundary_u),
							    IDATA(isDOFBoundary_v),
							    DDATA(n),
							    DDATA(bc_f),
							    DDATA(bc_fpu),
							    DDATA(bc_fpv),
							    DDATA(f),
							    DDATA(fpu),
							    DDATA(fpv),
							    DDATA(df_du),
							    DDATA(df_dv),
							    DDATA(dfpu_dp),
							    DDATA(dfpv_dp),
							    DDATA(advectiveFlux),
							    DDATA(advectiveFluxpu),
							    DDATA(advectiveFluxpv),
							    DDATA(dadvectiveFlux_du_left),
							    DDATA(dadvectiveFlux_dv_left),
							    DDATA(dadvectiveFluxpu_dp_left),
							    DDATA(dadvectiveFluxpv_dp_left));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokesP3D(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *bc_f,
    *bc_fpu,
    *bc_fpv,
    *bc_fpw,
    *f,
    *fpu,
    *fpv,
    *fpw,
    *df_du,
    *df_dv,
    *df_dw,
    *dfpu_dp,
    *dfpv_dp,
    *dfpw_dp,
    *advectiveFlux,
    *advectiveFluxpu,
    *advectiveFluxpv,
    *advectiveFluxpw,
    *dadvectiveFlux_du_left,
    *dadvectiveFlux_dv_left,
    *dadvectiveFlux_dw_left,
    *dadvectiveFluxpu_dp_left,
    *dadvectiveFluxpv_dp_left,
    *dadvectiveFluxpw_dp_left,
    *isDOFBoundary_p,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *isDOFBoundary_w;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_p,
                       &isDOFBoundary_u,
                       &isDOFBoundary_v,
                       &isDOFBoundary_w,
                       &n,
                       &bc_f,
                       &bc_fpu,
                       &bc_fpv,
                       &bc_fpw,
                       &f,
                       &fpu,
                       &fpv,
                       &fpw,
                       &df_du,
                       &df_dv,
                       &df_dw,
                       &dfpu_dp,
                       &dfpv_dp,
                       &dfpw_dp,
                       &advectiveFlux,
                       &advectiveFluxpu,
                       &advectiveFluxpv,
                       &advectiveFluxpw,
                       &dadvectiveFlux_du_left,
                       &dadvectiveFlux_dv_left,
                       &dadvectiveFlux_dw_left,
                       &dadvectiveFluxpu_dp_left,
                       &dadvectiveFluxpv_dp_left,
                       &dadvectiveFluxpw_dp_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFluxStokesP3D(SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(f)[1],
						      SHAPE(f)[2],
						      SHAPE(f)[3],
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(isDOFBoundary_p),
						      IDATA(isDOFBoundary_u),
						      IDATA(isDOFBoundary_v),
						      IDATA(isDOFBoundary_w),
						      DDATA(n),
						      DDATA(bc_f),
						      DDATA(bc_fpu),
						      DDATA(bc_fpv),
						      DDATA(bc_fpw),
						      DDATA(f),
						      DDATA(fpu),
						      DDATA(fpv),
						      DDATA(fpw),
						      DDATA(df_du),
						      DDATA(df_dv),
						      DDATA(df_dw),
						      DDATA(dfpu_dp),
						      DDATA(dfpv_dp),
						      DDATA(dfpw_dp),
						      DDATA(advectiveFlux),
						      DDATA(advectiveFluxpu),
						      DDATA(advectiveFluxpv),
						      DDATA(advectiveFluxpw),
						      DDATA(dadvectiveFlux_du_left),
						      DDATA(dadvectiveFlux_dv_left),
						      DDATA(dadvectiveFlux_dw_left),
						      DDATA(dadvectiveFluxpu_dp_left),
						      DDATA(dadvectiveFluxpv_dp_left),
						      DDATA(dadvectiveFluxpw_dp_left));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFluxStokesP3D(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(f)[1],
							    SHAPE(f)[2],
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary_p),
							    IDATA(isDOFBoundary_u),
							    IDATA(isDOFBoundary_v),
							    IDATA(isDOFBoundary_w),
							    DDATA(n),
							    DDATA(bc_f),
							    DDATA(bc_fpu),
							    DDATA(bc_fpv),
							    DDATA(bc_fpw),
							    DDATA(f),
							    DDATA(fpu),
							    DDATA(fpv),
							    DDATA(fpw),
							    DDATA(df_du),
							    DDATA(df_dv),
							    DDATA(df_dw),
							    DDATA(dfpu_dp),
							    DDATA(dfpv_dp),
							    DDATA(dfpw_dp),
							    DDATA(advectiveFlux),
							    DDATA(advectiveFluxpu),
							    DDATA(advectiveFluxpv),
							    DDATA(advectiveFluxpw),
							    DDATA(dadvectiveFlux_du_left),
							    DDATA(dadvectiveFlux_dv_left),
							    DDATA(dadvectiveFlux_dw_left),
							    DDATA(dadvectiveFluxpu_dp_left),
							    DDATA(dadvectiveFluxpv_dp_left),
							    DDATA(dadvectiveFluxpw_dp_left));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxNavierStokes2D(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *bc_p,
    *bc_f_mass,
    *bc_f_umom,
    *bc_f_vmom,
    *p,
    *dm_umom,
    *f_mass,
    *f_umom,
    *f_vmom,
    *df_mass_du,
    *df_mass_dv,
    *df_umom_dp,
    *df_umom_du,
    *df_umom_dv,
    *df_vmom_dp,
    *df_vmom_du,
    *df_vmom_dv,
    *advectiveFlux_mass,
    *advectiveFlux_umom,
    *advectiveFlux_vmom,
    *dadvectiveFlux_mass_dp,
    *dadvectiveFlux_mass_du,
    *dadvectiveFlux_mass_dv,
    *dadvectiveFlux_umom_dp,
    *dadvectiveFlux_umom_du,
    *dadvectiveFlux_umom_dv,
    *dadvectiveFlux_vmom_dp,
    *dadvectiveFlux_vmom_du,
    *dadvectiveFlux_vmom_dv,
    *isDOFBoundary_p,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *velocity;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_p,
                       &isDOFBoundary_u,
                       &isDOFBoundary_v,
                       &n,
                       &bc_p,
                       &bc_f_mass,
                       &bc_f_umom,
                       &bc_f_vmom,
                       &p,
                       &dm_umom,
                       &f_mass,
                       &f_umom,
                       &f_vmom,
                       &df_mass_du,
                       &df_mass_dv,
                       &df_umom_dp,
                       &df_umom_du,
                       &df_umom_dv,
                       &df_vmom_dp,
                       &df_vmom_du,
                       &df_vmom_dv,
                       &advectiveFlux_mass,
                       &advectiveFlux_umom,
                       &advectiveFlux_vmom,
                       &dadvectiveFlux_mass_dp,
                       &dadvectiveFlux_mass_du,
                       &dadvectiveFlux_mass_dv,
                       &dadvectiveFlux_umom_dp,
                       &dadvectiveFlux_umom_du,
                       &dadvectiveFlux_umom_dv,
                       &dadvectiveFlux_vmom_dp,
                       &dadvectiveFlux_vmom_du,
                       &dadvectiveFlux_vmom_dv,
		       &velocity))
    return NULL;
  if (ND(f_mass) > 3)
    {
      calculateExteriorNumericalAdvectiveFluxNavierStokes2D(SHAPE(exteriorElementBoundaries)[0],
                                                            SHAPE(f_mass)[1],
                                                            SHAPE(f_mass)[2],
                                                            SHAPE(f_mass)[3],
                                                            IDATA(exteriorElementBoundaries),
                                                            IDATA(elementBoundaryElements),
                                                            IDATA(elementBoundaryLocalElementBoundaries),
                                                            IDATA(isDOFBoundary_p),
                                                            IDATA(isDOFBoundary_u),
                                                            IDATA(isDOFBoundary_v),
                                                            DDATA(n),
                                                            DDATA(bc_p),
                                                            DDATA(bc_f_mass),
                                                            DDATA(bc_f_umom),
                                                            DDATA(bc_f_vmom),
                                                            DDATA(p),
                                                            DDATA(f_mass),
                                                            DDATA(f_umom),
                                                            DDATA(f_vmom),
                                                            DDATA(df_mass_du),
                                                            DDATA(df_mass_dv),
                                                            DDATA(df_umom_du),
                                                            DDATA(df_umom_dv),
                                                            DDATA(df_vmom_du),
                                                            DDATA(df_vmom_dv),
                                                            DDATA(advectiveFlux_mass),
                                                            DDATA(advectiveFlux_umom),
                                                            DDATA(advectiveFlux_vmom),
                                                            DDATA(dadvectiveFlux_mass_du),
                                                            DDATA(dadvectiveFlux_mass_dv),
                                                            DDATA(dadvectiveFlux_umom_dp),
                                                            DDATA(dadvectiveFlux_umom_du),
                                                            DDATA(dadvectiveFlux_umom_dv),
                                                            DDATA(dadvectiveFlux_vmom_dp),
                                                            DDATA(dadvectiveFlux_vmom_du),
                                                            DDATA(dadvectiveFlux_vmom_dv));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFluxNavierStokes2D(SHAPE(exteriorElementBoundaries)[0],
                                                                  SHAPE(f_mass)[1],
                                                                  SHAPE(f_mass)[2],
                                                                  IDATA(exteriorElementBoundaries),
                                                                  IDATA(elementBoundaryElements),
                                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                                  IDATA(isDOFBoundary_p),
                                                                  IDATA(isDOFBoundary_u),
                                                                  IDATA(isDOFBoundary_v),
                                                                  DDATA(n),
                                                                  DDATA(bc_p),
                                                                  DDATA(bc_f_mass),
                                                                  DDATA(bc_f_umom),
                                                                  DDATA(bc_f_vmom),
                                                                  DDATA(p),
                                                                  DDATA(dm_umom),
                                                                  DDATA(f_mass),
                                                                  DDATA(f_umom),
                                                                  DDATA(f_vmom),
                                                                  DDATA(df_mass_du),
                                                                  DDATA(df_mass_dv),
                                                                  DDATA(df_umom_dp),
                                                                  DDATA(df_umom_du),
                                                                  DDATA(df_umom_dv),
                                                                  DDATA(df_vmom_dp),
                                                                  DDATA(df_vmom_du),
                                                                  DDATA(df_vmom_dv),
                                                                  DDATA(advectiveFlux_mass),
                                                                  DDATA(advectiveFlux_umom),
                                                                  DDATA(advectiveFlux_vmom),
                                                                  DDATA(dadvectiveFlux_mass_dp),
                                                                  DDATA(dadvectiveFlux_mass_du),
                                                                  DDATA(dadvectiveFlux_mass_dv),
                                                                  DDATA(dadvectiveFlux_umom_dp),
                                                                  DDATA(dadvectiveFlux_umom_du),
                                                                  DDATA(dadvectiveFlux_umom_dv),
                                                                  DDATA(dadvectiveFlux_vmom_dp),
                                                                  DDATA(dadvectiveFlux_vmom_du),
                                                                  DDATA(dadvectiveFlux_vmom_dv),
								  DDATA(velocity));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxNavierStokes3D(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *bc_p,
    *bc_f_mass,
    *bc_f_umom,
    *bc_f_vmom,
    *bc_f_wmom,
    *p,
    *f_mass,
    *f_umom,
    *f_vmom,
    *f_wmom,
    *df_mass_du,
    *df_mass_dv,
    *df_mass_dw,
    *df_umom_dp,
    *df_umom_du,
    *df_umom_dv,
    *df_umom_dw,
    *df_vmom_dp,
    *df_vmom_du,
    *df_vmom_dv,
    *df_vmom_dw,
    *df_wmom_dp,
    *df_wmom_du,
    *df_wmom_dv,
    *df_wmom_dw,
    *advectiveFlux_mass,
    *advectiveFlux_umom,
    *advectiveFlux_vmom,
    *advectiveFlux_wmom,
    *dadvectiveFlux_mass_du,
    *dadvectiveFlux_mass_dv,
    *dadvectiveFlux_mass_dw,
    *dadvectiveFlux_umom_dp,
    *dadvectiveFlux_umom_du,
    *dadvectiveFlux_umom_dv,
    *dadvectiveFlux_umom_dw,
    *dadvectiveFlux_vmom_dp,
    *dadvectiveFlux_vmom_du,
    *dadvectiveFlux_vmom_dv,
    *dadvectiveFlux_vmom_dw,
    *dadvectiveFlux_wmom_dp,
    *dadvectiveFlux_wmom_du,
    *dadvectiveFlux_wmom_dv,
    *dadvectiveFlux_wmom_dw,
    *isDOFBoundary_p,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *isDOFBoundary_w,
    *velocity;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_p,
                       &isDOFBoundary_u,
                       &isDOFBoundary_v,
                       &isDOFBoundary_w,
                       &n,
                       &bc_p,
                       &bc_f_mass,
                       &bc_f_umom,
                       &bc_f_vmom,
                       &bc_f_wmom,
                       &p,
                       &f_mass,
                       &f_umom,
                       &f_vmom,
                       &f_wmom,
                       &df_mass_du,
                       &df_mass_dv,
                       &df_mass_dw,
                       &df_umom_dp,
                       &df_umom_du,
                       &df_umom_dv,
                       &df_umom_dw,
                       &df_vmom_dp,
                       &df_vmom_du,
                       &df_vmom_dv,
                       &df_vmom_dw,
                       &df_wmom_dp,
                       &df_wmom_du,
                       &df_wmom_dv,
                       &df_wmom_dw,
                       &advectiveFlux_mass,
                       &advectiveFlux_umom,
                       &advectiveFlux_vmom,
                       &advectiveFlux_wmom,
                       &dadvectiveFlux_mass_du,
                       &dadvectiveFlux_mass_dv,
                       &dadvectiveFlux_mass_dw,
                       &dadvectiveFlux_umom_dp,
                       &dadvectiveFlux_umom_du,
                       &dadvectiveFlux_umom_dv,
                       &dadvectiveFlux_umom_dw,
                       &dadvectiveFlux_vmom_dp,
                       &dadvectiveFlux_vmom_du,
                       &dadvectiveFlux_vmom_dv,
                       &dadvectiveFlux_vmom_dw,
                       &dadvectiveFlux_wmom_dp,
                       &dadvectiveFlux_wmom_du,
                       &dadvectiveFlux_wmom_dv,
                       &dadvectiveFlux_wmom_dw,
		       &velocity))
    return NULL;
  if (ND(f_mass) > 3)
    {
      exit(1);
/*       calculateExteriorNumericalAdvectiveFluxNavierStokes3D(SHAPE(exteriorElementBoundaries)[0], */
/*                                                             SHAPE(f_mass)[1], */
/*                                                             SHAPE(f_mass)[2], */
/*                                                             SHAPE(f_mass)[3], */
/*                                                             IDATA(exteriorElementBoundaries), */
/*                                                             IDATA(elementBoundaryElements), */
/*                                                             IDATA(elementBoundaryLocalElementBoundaries), */
/*                                                             IDATA(isDOFBoundary_p), */
/*                                                             IDATA(isDOFBoundary_u), */
/*                                                             IDATA(isDOFBoundary_v), */
/*                                                             DDATA(n), */
/*                                                             DDATA(bc_p), */
/*                                                             DDATA(bc_f_mass), */
/*                                                             DDATA(bc_f_umom), */
/*                                                             DDATA(bc_f_vmom), */
/*                                                             DDATA(p), */
/*                                                             DDATA(f_mass), */
/*                                                             DDATA(f_umom), */
/*                                                             DDATA(f_vmom), */
/*                                                             DDATA(df_mass_du), */
/*                                                             DDATA(df_mass_dv), */
/*                                                             DDATA(df_umom_du), */
/*                                                             DDATA(df_umom_dv), */
/*                                                             DDATA(df_vmom_du), */
/*                                                             DDATA(df_vmom_dv), */
/*                                                             DDATA(advectiveFlux_mass), */
/*                                                             DDATA(advectiveFlux_umom), */
/*                                                             DDATA(advectiveFlux_vmom), */
/*                                                             DDATA(dadvectiveFlux_mass_du), */
/*                                                             DDATA(dadvectiveFlux_mass_dv), */
/*                                                             DDATA(dadvectiveFlux_umom_dp), */
/*                                                             DDATA(dadvectiveFlux_umom_du), */
/*                                                             DDATA(dadvectiveFlux_umom_dv), */
/*                                                             DDATA(dadvectiveFlux_vmom_dp), */
/*                                                             DDATA(dadvectiveFlux_vmom_du), */
/*                                                             DDATA(dadvectiveFlux_vmom_dv)); */
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFluxNavierStokes3D(SHAPE(exteriorElementBoundaries)[0],
                                                                  SHAPE(f_mass)[1],
                                                                  SHAPE(f_mass)[2],
                                                                  IDATA(exteriorElementBoundaries),
                                                                  IDATA(elementBoundaryElements),
                                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                                  IDATA(isDOFBoundary_p),
                                                                  IDATA(isDOFBoundary_u),
                                                                  IDATA(isDOFBoundary_v),
                                                                  IDATA(isDOFBoundary_w),
                                                                  DDATA(n),
                                                                  DDATA(bc_p),
                                                                  DDATA(bc_f_mass),
                                                                  DDATA(bc_f_umom),
                                                                  DDATA(bc_f_vmom),
                                                                  DDATA(bc_f_wmom),
                                                                  DDATA(p),
                                                                  DDATA(f_mass),
                                                                  DDATA(f_umom),
                                                                  DDATA(f_vmom),
                                                                  DDATA(f_wmom),
                                                                  DDATA(df_mass_du),
                                                                  DDATA(df_mass_dv),
                                                                  DDATA(df_mass_dw),
                                                                  DDATA(df_umom_dp),
                                                                  DDATA(df_umom_du),
                                                                  DDATA(df_umom_dv),
                                                                  DDATA(df_umom_dw),
                                                                  DDATA(df_vmom_dp),
                                                                  DDATA(df_vmom_du),
                                                                  DDATA(df_vmom_dv),
                                                                  DDATA(df_vmom_dw),
                                                                  DDATA(df_wmom_dp),
                                                                  DDATA(df_wmom_du),
                                                                  DDATA(df_wmom_dv),
                                                                  DDATA(df_wmom_dw),
                                                                  DDATA(advectiveFlux_mass),
                                                                  DDATA(advectiveFlux_umom),
                                                                  DDATA(advectiveFlux_vmom),
                                                                  DDATA(advectiveFlux_wmom),
                                                                  DDATA(dadvectiveFlux_mass_du),
                                                                  DDATA(dadvectiveFlux_mass_dv),
                                                                  DDATA(dadvectiveFlux_mass_dw),
                                                                  DDATA(dadvectiveFlux_umom_dp),
                                                                  DDATA(dadvectiveFlux_umom_du),
                                                                  DDATA(dadvectiveFlux_umom_dv),
                                                                  DDATA(dadvectiveFlux_umom_dw),
                                                                  DDATA(dadvectiveFlux_vmom_dp),
                                                                  DDATA(dadvectiveFlux_vmom_du),
                                                                  DDATA(dadvectiveFlux_vmom_dv),
                                                                  DDATA(dadvectiveFlux_vmom_dw),
                                                                  DDATA(dadvectiveFlux_wmom_dp),
                                                                  DDATA(dadvectiveFlux_wmom_du),
                                                                  DDATA(dadvectiveFlux_wmom_dv),
                                                                  DDATA(dadvectiveFlux_wmom_dw),
								  DDATA(velocity));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokes2D(PyObject* self, 
                                                              PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *bc_p,
    *bc_f_mass,
    *p,
    *f_mass,
    *df_mass_du,
    *df_mass_dv,
    *advectiveFlux_mass,
    *advectiveFlux_umom,
    *advectiveFlux_vmom,
    *dadvectiveFlux_mass_du,
    *dadvectiveFlux_mass_dv,
    *dadvectiveFlux_umom_dp,
    *dadvectiveFlux_vmom_dp,
    *isDOFBoundary_p,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *velocity;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_p,
                       &isDOFBoundary_u,
                       &isDOFBoundary_v,
                       &n,
                       &bc_p,
                       &bc_f_mass,
                       &p,
                       &f_mass,
                       &df_mass_du,
                       &df_mass_dv,
                       &advectiveFlux_mass,
                       &advectiveFlux_umom,
                       &advectiveFlux_vmom,
                       &dadvectiveFlux_mass_du,
                       &dadvectiveFlux_mass_dv,
                       &dadvectiveFlux_umom_dp,
                       &dadvectiveFlux_vmom_dp,
		       &velocity))
    return NULL;
  calculateGlobalExteriorNumericalAdvectiveFluxStokes2D(SHAPE(exteriorElementBoundaries)[0],
                                                        SHAPE(f_mass)[1],
                                                        SHAPE(f_mass)[2],
                                                        IDATA(exteriorElementBoundaries),
                                                        IDATA(elementBoundaryElements),
                                                        IDATA(elementBoundaryLocalElementBoundaries),
                                                        IDATA(isDOFBoundary_p),
                                                        IDATA(isDOFBoundary_u),
                                                        IDATA(isDOFBoundary_v),
                                                        DDATA(n),
                                                        DDATA(bc_p),
                                                        DDATA(bc_f_mass),
                                                        DDATA(p),
                                                        DDATA(f_mass),
                                                        DDATA(df_mass_du),
                                                        DDATA(df_mass_dv),
                                                        DDATA(advectiveFlux_mass),
                                                        DDATA(advectiveFlux_umom),
                                                        DDATA(advectiveFlux_vmom),
                                                        DDATA(dadvectiveFlux_mass_du),
                                                        DDATA(dadvectiveFlux_mass_dv),
                                                        DDATA(dadvectiveFlux_umom_dp),
                                                        DDATA(dadvectiveFlux_vmom_dp),
                                                        DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokes3D(PyObject* self, 
                                                              PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *bc_p,
    *bc_f_mass,
    *bc_f_umom,
    *bc_f_vmom,
    *bc_f_wmom,
    *p,
    *f_mass,
    *f_umom,
    *f_vmom,
    *f_wmom,
    *df_mass_du,
    *df_mass_dv,
    *df_mass_dw,
    *advectiveFlux_mass,
    *advectiveFlux_umom,
    *advectiveFlux_vmom,
    *advectiveFlux_wmom,
    *dadvectiveFlux_mass_du,
    *dadvectiveFlux_mass_dv,
    *dadvectiveFlux_mass_dw,
    *dadvectiveFlux_umom_dp,
    *dadvectiveFlux_vmom_dp,
    *dadvectiveFlux_wmom_dp,
    *isDOFBoundary_p,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *isDOFBoundary_w,
    *velocity;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_p,
                       &isDOFBoundary_u,
                       &isDOFBoundary_v,
                       &isDOFBoundary_w,
                       &n,
                       &bc_p,
                       &bc_f_mass,
                       &p,
                       &f_mass,
                       &df_mass_du,
                       &df_mass_dv,
                       &df_mass_dw,
                       &advectiveFlux_mass,
                       &advectiveFlux_umom,
                       &advectiveFlux_vmom,
                       &advectiveFlux_wmom,
                       &dadvectiveFlux_mass_du,
                       &dadvectiveFlux_mass_dv,
                       &dadvectiveFlux_mass_dw,
                       &dadvectiveFlux_umom_dp,
                       &dadvectiveFlux_vmom_dp,
                       &dadvectiveFlux_wmom_dp,
		       &velocity))
    return NULL;
  calculateGlobalExteriorNumericalAdvectiveFluxStokes3D(SHAPE(exteriorElementBoundaries)[0],
                                                        SHAPE(f_mass)[1],
                                                        SHAPE(f_mass)[2],
                                                        IDATA(exteriorElementBoundaries),
                                                        IDATA(elementBoundaryElements),
                                                        IDATA(elementBoundaryLocalElementBoundaries),
                                                        IDATA(isDOFBoundary_p),
                                                        IDATA(isDOFBoundary_u),
                                                        IDATA(isDOFBoundary_v),
                                                        IDATA(isDOFBoundary_w),
                                                        DDATA(n),
                                                        DDATA(bc_p),
                                                        DDATA(bc_f_mass),
                                                        DDATA(p),
                                                        DDATA(f_mass),
                                                        DDATA(df_mass_du),
                                                        DDATA(df_mass_dv),
                                                        DDATA(df_mass_dw),
                                                        DDATA(advectiveFlux_mass),
                                                        DDATA(advectiveFlux_umom),
                                                        DDATA(advectiveFlux_vmom),
                                                        DDATA(advectiveFlux_wmom),
                                                        DDATA(dadvectiveFlux_mass_du),
                                                        DDATA(dadvectiveFlux_mass_dv),
                                                        DDATA(dadvectiveFlux_mass_dw),
                                                        DDATA(dadvectiveFlux_umom_dp),
                                                        DDATA(dadvectiveFlux_vmom_dp),
                                                        DDATA(dadvectiveFlux_wmom_dp),
                                                        DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorInflowNumericalAdvectiveFlux(PyObject* self, 
								  PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*f,*df,
    *advectiveFlux,*dadvectiveFlux_left,*inflowFlag,*inflowFlux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &inflowFlag,
		       &inflowFlux,
                       &n,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left))
    return NULL;
  calculateGlobalExteriorInflowNumericalAdvectiveFlux(SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(f)[1],
						      SHAPE(f)[2],
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(inflowFlag),
						      DDATA(inflowFlux),
						      DDATA(n),
						      DDATA(f),
						      DDATA(df),
						      DDATA(advectiveFlux),
						      DDATA(dadvectiveFlux_left));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalAdvectiveFlux_average(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*bc_u,*bc_f,*bc_df,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*inflowFlag,*isDOFBoundary;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &inflowFlag,
                       &n,
                       &bc_u,
                       &bc_f,
                       &bc_df,
                       &u,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left))
    return NULL;
  if (ND(f) > 3)
    {
      calculateExteriorNumericalAdvectiveFlux_average(SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(f)[1],
						      SHAPE(f)[2],
						      SHAPE(f)[3],
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(isDOFBoundary),
						      IDATA(inflowFlag),
						      DDATA(n),
						      DDATA(bc_u),
						      DDATA(bc_f),
						      DDATA(bc_df),
						      DDATA(u),
						      DDATA(f),
						      DDATA(df),
						      DDATA(advectiveFlux),
						      DDATA(dadvectiveFlux_left));
    }
  else
    {
      calculateGlobalExteriorNumericalAdvectiveFlux_average(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(f)[1],
							    SHAPE(f)[2],
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary),
							    IDATA(inflowFlag),
							    DDATA(n),
							    DDATA(bc_u),
							    DDATA(bc_f),
							    DDATA(bc_df),
							    DDATA(u),
							    DDATA(f),
							    DDATA(df),
							    DDATA(advectiveFlux),
							    DDATA(dadvectiveFlux_left));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalAdvectiveFluxJacobian(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *inflowFlag,
    *dadvectiveFlux_left,
    *v,
    *fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &inflowFlag,
                       &dadvectiveFlux_left,
		       &v,
                       &fluxJacobian))
    return NULL;
  if (ND(v) > 3)
    {
      assert(ND(v) == 4);
      updateExteriorNumericalAdvectiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(v)[1],
						   SHAPE(v)[2],
						   SHAPE(v)[3],
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(inflowFlag),
						   DDATA(dadvectiveFlux_left),
						   DDATA(v),
						   DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(v) == 3);
      updateGlobalExteriorNumericalAdvectiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(v)[1],
							 SHAPE(v)[2],
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(inflowFlag),
							 DDATA(dadvectiveFlux_left),
							 DDATA(v),
							 DDATA(fluxJacobian));
      
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalAdvectiveFluxJacobian_free(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *inflowFlag,
    *dadvectiveFlux_left,
    *v,
    *fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &inflowFlag,
                       &dadvectiveFlux_left,
		       &v,
                       &fluxJacobian))
    return NULL;
  if (ND(v) > 3)
    {
      assert(ND(v) == 4);
      updateExteriorNumericalAdvectiveFluxJacobian_free(SHAPE(exteriorElementBoundaries)[0],
                                                        SHAPE(v)[1],
                                                        SHAPE(v)[2],
                                                        SHAPE(v)[3],
                                                        IDATA(exteriorElementBoundaries),
                                                        IDATA(elementBoundaryElements),
                                                        IDATA(elementBoundaryLocalElementBoundaries),
                                                        IDATA(inflowFlag),
                                                        DDATA(dadvectiveFlux_left),
                                                        DDATA(v),
                                                        DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(v) == 3);
      updateGlobalExteriorNumericalAdvectiveFluxJacobian_free(SHAPE(exteriorElementBoundaries)[0],
                                                              SHAPE(v)[1],
                                                              SHAPE(v)[2],
                                                              IDATA(exteriorElementBoundaries),
                                                              IDATA(elementBoundaryElements),
                                                              IDATA(elementBoundaryLocalElementBoundaries),
                                                              IDATA(inflowFlag),
                                                              DDATA(dadvectiveFlux_left),
                                                              DDATA(v),
                                                              DDATA(fluxJacobian));
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalDiffusiveFlux(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*a,*grad_phi,*u,*penalty,*diffusiveFlux;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO|id",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux(scale_penalty,penalty_floor,
                                          SHAPE(interiorElementBoundaries)[0],
                                          SHAPE(grad_phi)[1],
                                          SHAPE(grad_phi)[2],
                                          SHAPE(grad_phi)[3],
                                          IDATA(interiorElementBoundaries),
                                          IDATA(elementBoundaryElements),
                                          IDATA(elementBoundaryLocalElementBoundaries),
                                          DDATA(n),
                                          DDATA(a),
                                          DDATA(grad_phi),
                                          DDATA(u),
                                          DDATA(penalty),
                                          DDATA(diffusiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}
/*
static PyObject*
cnumericalFluxCalculateInteriorNumericalDiffusiveFlux(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*a,*grad_phi,*u,*penalty,*diffusiveFlux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux(SHAPE(interiorElementBoundaries)[0],
                                          SHAPE(grad_phi)[1],
                                          SHAPE(grad_phi)[2],
                                          SHAPE(grad_phi)[3],
                                          IDATA(interiorElementBoundaries),
                                          IDATA(elementBoundaryElements),
                                          IDATA(elementBoundaryLocalElementBoundaries),
                                          DDATA(n),
                                          DDATA(a),
                                          DDATA(grad_phi),
                                          DDATA(u),
                                          DDATA(penalty),
                                          DDATA(diffusiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}
*/

static PyObject*
cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_sd(PyObject* self, 
							 PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*a,*grad_phi,*u,*penalty,*diffusiveFlux,*rowptr,*colind;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux_sd(scale_penalty,penalty_floor,
                                             SHAPE(interiorElementBoundaries)[0],
					     SHAPE(grad_phi)[1],
					     SHAPE(grad_phi)[2],
					     SHAPE(grad_phi)[3],
					     IDATA(rowptr),
					     IDATA(colind),
					     IDATA(interiorElementBoundaries),
					     IDATA(elementBoundaryElements),
					     IDATA(elementBoundaryLocalElementBoundaries),
					     DDATA(n),
					     DDATA(a),
					     DDATA(grad_phi),
					     DDATA(u),
					     DDATA(penalty),
					     DDATA(diffusiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}
/*
static PyObject*
cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_sd(PyObject* self, 
							 PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*a,*grad_phi,*u,*penalty,*diffusiveFlux,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux_sd(SHAPE(interiorElementBoundaries)[0],
					     SHAPE(grad_phi)[1],
					     SHAPE(grad_phi)[2],
					     SHAPE(grad_phi)[3],
					     IDATA(rowptr),
					     IDATA(colind),
					     IDATA(interiorElementBoundaries),
					     IDATA(elementBoundaryElements),
					     IDATA(elementBoundaryLocalElementBoundaries),
					     DDATA(n),
					     DDATA(a),
					     DDATA(grad_phi),
					     DDATA(u),
					     DDATA(penalty),
					     DDATA(diffusiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}
*/

static PyObject*
cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO|id",
                       &l2g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian(scale_penalty,penalty_floor,
                                               SHAPE(interiorElementBoundaries)[0],
                                               SHAPE(grad_v)[1],
                                               SHAPE(grad_v)[2],
                                               SHAPE(grad_v)[3],
                                               SHAPE(grad_v)[4],
                                               IDATA(l2g),
                                               IDATA(interiorElementBoundaries),
                                               IDATA(elementBoundaryElements),
                                               IDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(a),
                                               DDATA(da),
                                               DDATA(grad_phi),
                                               DDATA(dphi),
                                               DDATA(v),
                                               DDATA(grad_v),
                                               DDATA(penalty),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}
/*
static PyObject*
cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
                       &l2g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian(SHAPE(interiorElementBoundaries)[0],
                                               SHAPE(grad_v)[1],
                                               SHAPE(grad_v)[2],
                                               SHAPE(grad_v)[3],
                                               SHAPE(grad_v)[4],
                                               IDATA(l2g),
                                               IDATA(interiorElementBoundaries),
                                               IDATA(elementBoundaryElements),
                                               IDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(a),
                                               DDATA(da),
                                               DDATA(grad_phi),
                                               DDATA(dphi),
                                               DDATA(v),
                                               DDATA(grad_v),
                                               DDATA(penalty),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}
*/

static PyObject*
cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_sd(PyObject* self, 
							      PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian,*rowptr,*colind;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
                       &l2g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian_sd(scale_penalty,penalty_floor,
                                                  SHAPE(interiorElementBoundaries)[0],
						  SHAPE(grad_v)[1],
						  SHAPE(grad_v)[2],
						  SHAPE(grad_v)[3],
						  SHAPE(grad_v)[4],
						  IDATA(rowptr),
						  IDATA(colind),
						  IDATA(l2g),
						  IDATA(interiorElementBoundaries),
						  IDATA(elementBoundaryElements),
						  IDATA(elementBoundaryLocalElementBoundaries),
						  DDATA(n),
						  DDATA(a),
						  DDATA(da),
						  DDATA(grad_phi),
						  DDATA(dphi),
						  DDATA(v),
						  DDATA(grad_v),
						  DDATA(penalty),
						  DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}
/*
static PyObject*
cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_sd(PyObject* self, 
							      PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &l2g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian_sd(SHAPE(interiorElementBoundaries)[0],
						  SHAPE(grad_v)[1],
						  SHAPE(grad_v)[2],
						  SHAPE(grad_v)[3],
						  SHAPE(grad_v)[4],
						  IDATA(rowptr),
						  IDATA(colind),
						  IDATA(l2g),
						  IDATA(interiorElementBoundaries),
						  IDATA(elementBoundaryElements),
						  IDATA(elementBoundaryLocalElementBoundaries),
						  DDATA(n),
						  DDATA(a),
						  DDATA(da),
						  DDATA(grad_phi),
						  DDATA(dphi),
						  DDATA(v),
						  DDATA(grad_v),
						  DDATA(penalty),
						  DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}
*/
static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO|id",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      calculateExteriorNumericalDiffusiveFlux(scale_penalty,penalty_floor,
                                              SHAPE(exteriorElementBoundaries)[0],
					      SHAPE(grad_phi)[1],
					      SHAPE(grad_phi)[2],
					      SHAPE(grad_phi)[3],
					      IDATA(exteriorElementBoundaries),
					      IDATA(elementBoundaryElements),
					      IDATA(elementBoundaryLocalElementBoundaries),
					      IDATA(isDOFBoundary),
					      DDATA(n),
					      DDATA(bc_a),
					      DDATA(bc_grad_phi),
					      DDATA(bc_u),
					      DDATA(a),
					      DDATA(grad_phi),
					      DDATA(u),
					      DDATA(penalty),
					      DDATA(diffusiveFlux));
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux(scale_penalty,penalty_floor,
                                                    SHAPE(exteriorElementBoundaries)[0],
						    SHAPE(grad_phi)[1],
						    SHAPE(grad_phi)[2],
						    IDATA(exteriorElementBoundaries),
						    IDATA(elementBoundaryElements),
						    IDATA(elementBoundaryLocalElementBoundaries),
						    IDATA(isDOFBoundary),
						    DDATA(n),
						    DDATA(bc_a),
						    DDATA(bc_grad_phi),
						    DDATA(bc_u),
						    DDATA(a),
						    DDATA(grad_phi),
						    DDATA(u),
						    DDATA(penalty),
						    DDATA(diffusiveFlux));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

/*
static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      calculateExteriorNumericalDiffusiveFlux(SHAPE(exteriorElementBoundaries)[0],
					      SHAPE(grad_phi)[1],
					      SHAPE(grad_phi)[2],
					      SHAPE(grad_phi)[3],
					      IDATA(exteriorElementBoundaries),
					      IDATA(elementBoundaryElements),
					      IDATA(elementBoundaryLocalElementBoundaries),
					      IDATA(isDOFBoundary),
					      DDATA(n),
					      DDATA(bc_a),
					      DDATA(bc_grad_phi),
					      DDATA(bc_u),
					      DDATA(a),
					      DDATA(grad_phi),
					      DDATA(u),
					      DDATA(penalty),
					      DDATA(diffusiveFlux));
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux(SHAPE(exteriorElementBoundaries)[0],
						    SHAPE(grad_phi)[1],
						    SHAPE(grad_phi)[2],
						    IDATA(exteriorElementBoundaries),
						    IDATA(elementBoundaryElements),
						    IDATA(elementBoundaryLocalElementBoundaries),
						    IDATA(isDOFBoundary),
						    DDATA(n),
						    DDATA(bc_a),
						    DDATA(bc_grad_phi),
						    DDATA(bc_u),
						    DDATA(a),
						    DDATA(grad_phi),
						    DDATA(u),
						    DDATA(penalty),
						    DDATA(diffusiveFlux));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}
*/

static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_sd(PyObject* self, 
							 PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux,*rowptr,*colind;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      calculateExteriorNumericalDiffusiveFlux_sd(scale_penalty,penalty_floor,
                                                 SHAPE(exteriorElementBoundaries)[0],
						 SHAPE(grad_phi)[1],
						 SHAPE(grad_phi)[2],
						 SHAPE(grad_phi)[3],
						 IDATA(rowptr),
						 IDATA(colind),
						 IDATA(exteriorElementBoundaries),
						 IDATA(elementBoundaryElements),
						 IDATA(elementBoundaryLocalElementBoundaries),
						 IDATA(isDOFBoundary),
						 DDATA(n),
						 DDATA(bc_a),
						 DDATA(bc_grad_phi),
						 DDATA(bc_u),
						 DDATA(a),
						 DDATA(grad_phi),
						 DDATA(u),
						 DDATA(penalty),
						 DDATA(diffusiveFlux));
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux_sd(scale_penalty,penalty_floor,
                                                       SHAPE(exteriorElementBoundaries)[0],
						       SHAPE(grad_phi)[1],
						       SHAPE(grad_phi)[2],
						       IDATA(rowptr),
						       IDATA(colind),
						       IDATA(exteriorElementBoundaries),
						       IDATA(elementBoundaryElements),
						       IDATA(elementBoundaryLocalElementBoundaries),
						       IDATA(isDOFBoundary),
						       DDATA(n),
						       DDATA(bc_a),
						       DDATA(bc_grad_phi),
						       DDATA(bc_u),
						       DDATA(a),
						       DDATA(grad_phi),
						       DDATA(u),
						       DDATA(penalty),
						       DDATA(diffusiveFlux));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFluxWithUpwinding_sd(PyObject* self, 
								      PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*fluxBoundaryFlag,
    *n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux,*rowptr,*colind;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
		       &fluxBoundaryFlag,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      assert(0); /*not implemented*/
      /* calculateExteriorNumericalDiffusiveFluxWithUpwinding_sd(scale_penalty,penalty_floor, */
      /* 							      SHAPE(exteriorElementBoundaries)[0], */
      /* 							      SHAPE(grad_phi)[1], */
      /* 							      SHAPE(grad_phi)[2], */
      /* 							      SHAPE(grad_phi)[3], */
      /* 							      IDATA(rowptr), */
      /* 							      IDATA(colind), */
      /* 							      IDATA(exteriorElementBoundaries), */
      /* 							      IDATA(elementBoundaryElements), */
      /* 							      IDATA(elementBoundaryLocalElementBoundaries), */
      /* 							      IDATA(isDOFBoundary), */
      /* 							      IDATA(fluxBoundaryFlag), */
      /* 							      DDATA(n), */
      /* 							      DDATA(bc_a), */
      /* 							      DDATA(bc_grad_phi), */
      /* 							      DDATA(bc_u), */
      /* 							      DDATA(a), */
      /* 							      DDATA(grad_phi), */
      /* 							      DDATA(u), */
      /* 							      DDATA(penalty), */
      /* 							      DDATA(diffusiveFlux)); */
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux_sd(scale_penalty,penalty_floor,
                                                       SHAPE(exteriorElementBoundaries)[0],
						       SHAPE(grad_phi)[1],
						       SHAPE(grad_phi)[2],
						       IDATA(rowptr),
						       IDATA(colind),
						       IDATA(exteriorElementBoundaries),
						       IDATA(elementBoundaryElements),
						       IDATA(elementBoundaryLocalElementBoundaries),
						       IDATA(isDOFBoundary),
						       IDATA(fluxBoundaryFlag),
						       DDATA(n),
						       DDATA(bc_a),
						       DDATA(bc_grad_phi),
						       DDATA(bc_u),
						       DDATA(a),
						       DDATA(grad_phi),
						       DDATA(u),
						       DDATA(penalty),
						       DDATA(diffusiveFlux));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

/*
static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_sd(PyObject* self, 
							 PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      calculateExteriorNumericalDiffusiveFlux_sd(SHAPE(exteriorElementBoundaries)[0],
						 SHAPE(grad_phi)[1],
						 SHAPE(grad_phi)[2],
						 SHAPE(grad_phi)[3],
						 IDATA(rowptr),
						 IDATA(colind),
						 IDATA(exteriorElementBoundaries),
						 IDATA(elementBoundaryElements),
						 IDATA(elementBoundaryLocalElementBoundaries),
						 IDATA(isDOFBoundary),
						 DDATA(n),
						 DDATA(bc_a),
						 DDATA(bc_grad_phi),
						 DDATA(bc_u),
						 DDATA(a),
						 DDATA(grad_phi),
						 DDATA(u),
						 DDATA(penalty),
						 DDATA(diffusiveFlux));
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux_sd(SHAPE(exteriorElementBoundaries)[0],
						       SHAPE(grad_phi)[1],
						       SHAPE(grad_phi)[2],
						       IDATA(rowptr),
						       IDATA(colind),
						       IDATA(exteriorElementBoundaries),
						       IDATA(elementBoundaryElements),
						       IDATA(elementBoundaryLocalElementBoundaries),
						       IDATA(isDOFBoundary),
						       DDATA(n),
						       DDATA(bc_a),
						       DDATA(bc_grad_phi),
						       DDATA(bc_u),
						       DDATA(a),
						       DDATA(grad_phi),
						       DDATA(u),
						       DDATA(penalty),
						       DDATA(diffusiveFlux));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}
*/
static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_free(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      calculateExteriorNumericalDiffusiveFlux_free(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(grad_phi)[1],
						   SHAPE(grad_phi)[2],
						   SHAPE(grad_phi)[3],
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary),
						   DDATA(n),
						   DDATA(bc_a),
						   DDATA(bc_grad_phi),
						   DDATA(bc_u),
						   DDATA(a),
						   DDATA(grad_phi),
						   DDATA(u),
						   DDATA(penalty),
						   DDATA(diffusiveFlux));
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux_free(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(grad_phi)[1],
							 SHAPE(grad_phi)[2],
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(isDOFBoundary),
							 DDATA(n),
							 DDATA(bc_a),
							 DDATA(bc_grad_phi),
							 DDATA(bc_u),
							 DDATA(a),
							 DDATA(grad_phi),
							 DDATA(u),
							 DDATA(penalty),
							 DDATA(diffusiveFlux));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_free_sd(PyObject* self, 
							      PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_u,*a,*grad_phi,*u,*penalty,*diffusiveFlux,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_u,
                       &a,
                       &grad_phi,
                       &u,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == 4);
      calculateExteriorNumericalDiffusiveFlux_free_sd(SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(grad_phi)[1],
						      SHAPE(grad_phi)[2],
						      SHAPE(grad_phi)[3],
						      IDATA(rowptr),
						      IDATA(colind),
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(isDOFBoundary),
						      DDATA(n),
						      DDATA(bc_a),
						      DDATA(bc_grad_phi),
						      DDATA(bc_u),
						      DDATA(a),
						      DDATA(grad_phi),
						      DDATA(u),
						      DDATA(penalty),
						      DDATA(diffusiveFlux));
    }
  else
    {
      assert(ND(grad_phi) == 3);
      calculateGlobalExteriorNumericalDiffusiveFlux_free_sd(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(grad_phi)[1],
							    SHAPE(grad_phi)[2],
							    IDATA(rowptr),
							    IDATA(colind),
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary),
							    DDATA(n),
							    DDATA(bc_a),
							    DDATA(bc_grad_phi),
							    DDATA(bc_u),
							    DDATA(a),
							    DDATA(grad_phi),
							    DDATA(u),
							    DDATA(penalty),
							    DDATA(diffusiveFlux));

    }
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO|id",
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      updateExteriorNumericalDiffusiveFluxJacobian(scale_penalty,penalty_floor,
                                                   SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(grad_v)[1],
						   SHAPE(grad_v)[2],
						   SHAPE(grad_v)[3],
						   SHAPE(grad_v)[4],
						   IDATA(l2g),
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary),
						   DDATA(n),
						   DDATA(a),
						   DDATA(da),
						   DDATA(grad_phi),
						   DDATA(dphi),
						   DDATA(v),
						   DDATA(grad_v),
						   DDATA(penalty),
						   DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxJacobian(scale_penalty,penalty_floor,
                                                         SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(grad_v)[1],
							 SHAPE(grad_v)[2],
							 SHAPE(grad_v)[3],
							 IDATA(l2g),
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(isDOFBoundary),
							 DDATA(n),
							 DDATA(a),
							 DDATA(da),
							 DDATA(grad_phi),
							 DDATA(dphi),
							 DDATA(v),
							 DDATA(grad_v),
							 DDATA(penalty),
							 DDATA(fluxJacobian));
    }
  Py_INCREF(Py_None);
  return Py_None;
}

/*static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      updateExteriorNumericalDiffusiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(grad_v)[1],
						   SHAPE(grad_v)[2],
						   SHAPE(grad_v)[3],
						   SHAPE(grad_v)[4],
						   IDATA(l2g),
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary),
						   DDATA(n),
						   DDATA(a),
						   DDATA(da),
						   DDATA(grad_phi),
						   DDATA(dphi),
						   DDATA(v),
						   DDATA(grad_v),
						   DDATA(penalty),
						   DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(grad_v)[1],
							 SHAPE(grad_v)[2],
							 SHAPE(grad_v)[3],
							 IDATA(l2g),
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(isDOFBoundary),
							 DDATA(n),
							 DDATA(a),
							 DDATA(da),
							 DDATA(grad_phi),
							 DDATA(dphi),
							 DDATA(v),
							 DDATA(grad_v),
							 DDATA(penalty),
							 DDATA(fluxJacobian));
    }
  Py_INCREF(Py_None);
  return Py_None;
}
*/

static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_sd(PyObject* self, 
							      PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian,*rowptr,*colind;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      updateExteriorNumericalDiffusiveFluxJacobian_sd(scale_penalty,penalty_floor,
                                                      SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(grad_v)[1],
						      SHAPE(grad_v)[2],
						      SHAPE(grad_v)[3],
						      SHAPE(grad_v)[4],
						      IDATA(rowptr),
						      IDATA(colind),
						      IDATA(l2g),
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(isDOFBoundary),
						      DDATA(n),
						      DDATA(a),
						      DDATA(da),
						      DDATA(grad_phi),
						      DDATA(dphi),
						      DDATA(v),
						      DDATA(grad_v),
						      DDATA(penalty),
						      DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxJacobian_sd(scale_penalty,penalty_floor,
                                                            SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(grad_v)[1],
							    SHAPE(grad_v)[2],
							    SHAPE(grad_v)[3],
							    IDATA(rowptr),
							    IDATA(colind),
							    IDATA(l2g),
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary),
							    DDATA(n),
							    DDATA(a),
							    DDATA(da),
							    DDATA(grad_phi),
							    DDATA(dphi),
							    DDATA(v),
							    DDATA(grad_v),
							    DDATA(penalty),
							    DDATA(fluxJacobian));
    }
  Py_INCREF(Py_None);
  return Py_None;
}

/* 
static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_sd(PyObject* self, 
							      PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      updateExteriorNumericalDiffusiveFluxJacobian_sd(SHAPE(exteriorElementBoundaries)[0],
						      SHAPE(grad_v)[1],
						      SHAPE(grad_v)[2],
						      SHAPE(grad_v)[3],
						      SHAPE(grad_v)[4],
						      IDATA(rowptr),
						      IDATA(colind),
						      IDATA(l2g),
						      IDATA(exteriorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      IDATA(isDOFBoundary),
						      DDATA(n),
						      DDATA(a),
						      DDATA(da),
						      DDATA(grad_phi),
						      DDATA(dphi),
						      DDATA(v),
						      DDATA(grad_v),
						      DDATA(penalty),
						      DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxJacobian_sd(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(grad_v)[1],
							    SHAPE(grad_v)[2],
							    SHAPE(grad_v)[3],
							    IDATA(rowptr),
							    IDATA(colind),
							    IDATA(l2g),
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary),
							    DDATA(n),
							    DDATA(a),
							    DDATA(da),
							    DDATA(grad_phi),
							    DDATA(dphi),
							    DDATA(v),
							    DDATA(grad_v),
							    DDATA(penalty),
							    DDATA(fluxJacobian));
    }
  Py_INCREF(Py_None);
  return Py_None;
}
*/

static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_free(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      updateExteriorNumericalDiffusiveFluxJacobian_free(SHAPE(exteriorElementBoundaries)[0],
							SHAPE(grad_v)[1],
							SHAPE(grad_v)[2],
							SHAPE(grad_v)[3],
							SHAPE(grad_v)[4],
							IDATA(l2g),
							IDATA(exteriorElementBoundaries),
							IDATA(elementBoundaryElements),
							IDATA(elementBoundaryLocalElementBoundaries),
							IDATA(isDOFBoundary),
							DDATA(n),
							DDATA(a),
							DDATA(da),
							DDATA(grad_phi),
							DDATA(dphi),
							DDATA(v),
							DDATA(grad_v),
							DDATA(penalty),
							DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxJacobian_free(SHAPE(exteriorElementBoundaries)[0],
							      SHAPE(grad_v)[1],
							      SHAPE(grad_v)[2],
							      SHAPE(grad_v)[3],
							      IDATA(l2g),
							      IDATA(exteriorElementBoundaries),
							      IDATA(elementBoundaryElements),
							      IDATA(elementBoundaryLocalElementBoundaries),
							      IDATA(isDOFBoundary),
							      DDATA(n),
							      DDATA(a),
							      DDATA(da),
							      DDATA(grad_phi),
							      DDATA(dphi),
							      DDATA(v),
							      DDATA(grad_v),
							      DDATA(penalty),
							      DDATA(fluxJacobian));
      
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_free_sd(PyObject* self, 
								   PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      updateExteriorNumericalDiffusiveFluxJacobian_free_sd(SHAPE(exteriorElementBoundaries)[0],
							   SHAPE(grad_v)[1],
							   SHAPE(grad_v)[2],
							   SHAPE(grad_v)[3],
							   SHAPE(grad_v)[4],
							   IDATA(rowptr),
							   IDATA(colind),
							   IDATA(l2g),
							   IDATA(exteriorElementBoundaries),
							   IDATA(elementBoundaryElements),
							   IDATA(elementBoundaryLocalElementBoundaries),
							   IDATA(isDOFBoundary),
							   DDATA(n),
							   DDATA(a),
							   DDATA(da),
							   DDATA(grad_phi),
							   DDATA(dphi),
							   DDATA(v),
							   DDATA(grad_v),
							   DDATA(penalty),
							   DDATA(fluxJacobian));
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxJacobian_free_sd(SHAPE(exteriorElementBoundaries)[0],
								 SHAPE(grad_v)[1],
								 SHAPE(grad_v)[2],
								 SHAPE(grad_v)[3],
								 IDATA(rowptr),
								 IDATA(colind),
								 IDATA(l2g),
								 IDATA(exteriorElementBoundaries),
								 IDATA(elementBoundaryElements),
								 IDATA(elementBoundaryLocalElementBoundaries),
								 IDATA(isDOFBoundary),
								 DDATA(n),
								 DDATA(a),
								 DDATA(da),
								 DDATA(grad_phi),
								 DDATA(dphi),
								 DDATA(v),
								 DDATA(grad_v),
								 DDATA(penalty),
								 DDATA(fluxJacobian));
      
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxWithUpwindingJacobian_sd(PyObject* self, 
									   PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*l2g,*isDOFBoundary,*fluxBoundaryFlag,
    *n,*a,*da,*grad_phi,*dphi,*v,*grad_v,*penalty,*fluxJacobian,*rowptr,*colind;
  int scale_penalty = 0;
  double penalty_floor = 0.0;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOO|id",
		       &rowptr,
		       &colind,
                       &l2g,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &fluxBoundaryFlag,
                       &n,
                       &a,
                       &da,
                       &grad_phi,
                       &dphi,
                       &v,
                       &grad_v,
                       &penalty,
                       &fluxJacobian,
                       &scale_penalty,
                       &penalty_floor))
    return NULL;
  if (ND(grad_v) > 4)
    {
      assert(ND(grad_v) == 5);
      assert(0); /*not implemented*/
      /* updateExteriorNumericalDiffusiveFluxWithUpwindingJacobian_sd(scale_penalty,penalty_floor, */
      /* 								   SHAPE(exteriorElementBoundaries)[0], */
      /* 								   SHAPE(grad_v)[1], */
      /* 								   SHAPE(grad_v)[2], */
      /* 								   SHAPE(grad_v)[3], */
      /* 								   SHAPE(grad_v)[4], */
      /* 								   IDATA(rowptr), */
      /* 								   IDATA(colind), */
      /* 								   IDATA(l2g), */
      /* 								   IDATA(exteriorElementBoundaries), */
      /* 								   IDATA(elementBoundaryElements), */
      /* 								   IDATA(elementBoundaryLocalElementBoundaries), */
      /* 								   IDATA(isDOFBoundary), */
      /* 								   DDATA(n), */
      /* 								   DDATA(a), */
      /* 								   DDATA(da), */
      /* 								   DDATA(grad_phi), */
      /* 								   DDATA(dphi), */
      /* 								   DDATA(v), */
      /* 								   DDATA(grad_v), */
      /* 								   DDATA(penalty), */
      /* 								   DDATA(fluxJacobian)); */
    }
  else
    {
      assert(ND(grad_v) == 4);
      updateGlobalExteriorNumericalDiffusiveFluxWithUpwindingJacobian_sd(scale_penalty,penalty_floor,
									 SHAPE(exteriorElementBoundaries)[0],
									 SHAPE(grad_v)[1],
									 SHAPE(grad_v)[2],
									 SHAPE(grad_v)[3],
									 IDATA(rowptr),
									 IDATA(colind),
									 IDATA(l2g),
									 IDATA(exteriorElementBoundaries),
									 IDATA(elementBoundaryElements),
									 IDATA(elementBoundaryLocalElementBoundaries),
									 IDATA(isDOFBoundary),
									 DDATA(n),
									 DDATA(a),
									 DDATA(da),
									 DDATA(grad_phi),
									 DDATA(dphi),
									 DDATA(v),
									 DDATA(grad_v),
									 DDATA(penalty),
									 DDATA(fluxJacobian));
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalAdvectiveFlux(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*dadvectiveFlux_right;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left,
                       &dadvectiveFlux_right))
    return NULL;
  calculateInteriorNumericalAdvectiveFlux(SHAPE(interiorElementBoundaries)[0],
                                          SHAPE(f)[1],
                                          SHAPE(f)[2],
                                          SHAPE(f)[3],
                                          IDATA(interiorElementBoundaries),
                                          IDATA(elementBoundaryElements),
                                          IDATA(elementBoundaryLocalElementBoundaries),
                                          DDATA(n),
                                          DDATA(u),
                                          DDATA(f),
                                          DDATA(df),
                                          DDATA(advectiveFlux),
                                          DDATA(dadvectiveFlux_left),
                                          DDATA(dadvectiveFlux_right));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateInteriorNumericalAdvectiveFluxJacobian(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*dflux_left,*dflux_right,*v,*fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &dflux_left,
                       &dflux_right,
		       &v,
                       &fluxJacobian))
    return NULL;
  updateInteriorNumericalAdvectiveFluxJacobian(SHAPE(interiorElementBoundaries)[0],
                                               SHAPE(v)[1],
                                               SHAPE(v)[2],
                                               SHAPE(v)[3],
                                               IDATA(interiorElementBoundaries),
                                               IDATA(elementBoundaryElements),
                                               IDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(dflux_left),
                                               DDATA(dflux_right),
					       DDATA(v),
                                               DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxUpdateInteriorTwoSidedNumericalFluxJacobian(PyObject* self, 
							  PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*dflux_left,*dflux_right,*v,*fluxJacobian_2sided;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &dflux_left,
                       &dflux_right,
		       &v,
                       &fluxJacobian_2sided))
    return NULL;
  updateInteriorTwoSidedNumericalFluxJacobian(SHAPE(interiorElementBoundaries)[0],
					      SHAPE(v)[1],
					      SHAPE(v)[2],
					      SHAPE(v)[3],
					      IDATA(interiorElementBoundaries),
					      IDATA(elementBoundaryElements),
					      IDATA(elementBoundaryLocalElementBoundaries),
					      DDATA(dflux_left),
					      DDATA(dflux_right),
					      DDATA(v),
					      DDATA(fluxJacobian_2sided));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalAdvectiveFlux_average(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*u,*f,*df,*advectiveFlux,*dadvectiveFlux_left,*dadvectiveFlux_right;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &f,
                       &df,
                       &advectiveFlux,
                       &dadvectiveFlux_left,
                       &dadvectiveFlux_right))
    return NULL;
  calculateInteriorNumericalAdvectiveFlux_average(SHAPE(interiorElementBoundaries)[0],
                                                  SHAPE(f)[1],
                                                  SHAPE(f)[2],
                                                  SHAPE(f)[3],
                                                  IDATA(interiorElementBoundaries),
                                                  IDATA(elementBoundaryElements),
                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                  DDATA(n),
                                                  DDATA(u),
                                                  DDATA(f),
                                                  DDATA(df),
                                                  DDATA(advectiveFlux),
                                                  DDATA(dadvectiveFlux_left),
                                                  DDATA(dadvectiveFlux_right));
  Py_INCREF(Py_None);
  return Py_None;
}





static PyObject*
cnumericalFluxSetInflowFlux(PyObject* self, 
                           PyObject* args)
{
  PyObject *exteriorElementBoundaries,*flux,*inflowFlux;
  if(!PyArg_ParseTuple(args,"OOO",
                       &exteriorElementBoundaries,
                       &flux,
                       &inflowFlux))
    return NULL;
  setInflowFlux(SHAPE(exteriorElementBoundaries)[0],
                SHAPE(flux)[1],
                IDATA(exteriorElementBoundaries),
                DDATA(flux),
                DDATA(inflowFlux));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_LDG_upwind(PyObject* self,
                                                                PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *u,
    *a,
    *phi,
    *V,
    *penalty,
    *flux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &a,
                       &phi,
                       &V,
                       &penalty,
                       &flux))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux_LDG_upwind(SHAPE(interiorElementBoundaries)[0],
                                                     SHAPE(n)[1],
                                                     SHAPE(n)[2],
                                                     SHAPE(n)[3],
                                                     IDATA(interiorElementBoundaries),
                                                     IDATA(elementBoundaryElements),
                                                     IDATA(elementBoundaryLocalElementBoundaries),
                                                     DDATA(n),
                                                     DDATA(u),
                                                     DDATA(a),
                                                     DDATA(phi),
                                                     DDATA(V),
                                                     DDATA(penalty),
                                                     DDATA(flux));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_LDG_upwind_sd(PyObject* self,
								    PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *u,
    *a,
    *phi,
    *V,
    *penalty,
    *flux,
    *rowptr,
    *colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &a,
                       &phi,
                       &V,
                       &penalty,
                       &flux))
    return NULL;
  calculateInteriorNumericalDiffusiveFlux_LDG_upwind_sd(SHAPE(interiorElementBoundaries)[0],
							SHAPE(n)[1],
							SHAPE(n)[2],
							SHAPE(n)[3],
							IDATA(rowptr),
							IDATA(colind),
							IDATA(interiorElementBoundaries),
							IDATA(elementBoundaryElements),
							IDATA(elementBoundaryLocalElementBoundaries),
							DDATA(n),
							DDATA(u),
							DDATA(a),
							DDATA(phi),
							DDATA(V),
							DDATA(penalty),
							DDATA(flux));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind(PyObject* self,
                                                                        PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *a,
    *da,
    *dphi,
    *V,
    *dV,
    *dV_eb,
    *v,
    *penalty,
    *fluxJacobian,
    *fluxJacobian_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &dphi,
                       &V,
                       &dV,
                       &dV_eb,
                       &v,
                       &penalty,
                       &fluxJacobian,
                       &fluxJacobian_eb))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind(SHAPE(interiorElementBoundaries)[0],
                                                          SHAPE(v)[1],
                                                          SHAPE(v)[2],
                                                          SHAPE(v)[3],
                                                          SHAPE(n)[3],
                                                          IDATA(interiorElementBoundaries),
                                                          IDATA(elementBoundaryElements),
                                                          IDATA(elementBoundaryLocalElementBoundaries),
                                                          DDATA(n),
                                                          DDATA(a),
                                                          DDATA(da),
                                                          DDATA(dphi),
                                                          DDATA(V),
                                                          DDATA(dV),
                                                          DDATA(dV_eb),
                                                          DDATA(v),
                                                          DDATA(penalty),
                                                          DDATA(fluxJacobian),
                                                          DDATA(fluxJacobian_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(PyObject* self,
									 PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *a,
    *da,
    *dphi,
    *V,
    *dV,
    *dV_eb,
    *v,
    *penalty,
    *fluxJacobian,
    *fluxJacobian_eb,
    *rowptr,
    *colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &dphi,
                       &V,
                       &dV,
                       &dV_eb,
                       &v,
                       &penalty,
                       &fluxJacobian,
                       &fluxJacobian_eb))
    return NULL;
  updateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(SHAPE(interiorElementBoundaries)[0],
							     SHAPE(v)[1],
							     SHAPE(v)[2],
							     SHAPE(v)[3],
							     SHAPE(n)[3],
							     IDATA(rowptr),
							     IDATA(colind),
							     IDATA(interiorElementBoundaries),
							     IDATA(elementBoundaryElements),
							     IDATA(elementBoundaryLocalElementBoundaries),
							     DDATA(n),
							     DDATA(a),
							     DDATA(da),
							     DDATA(dphi),
							     DDATA(V),
							     DDATA(dV),
							     DDATA(dV_eb),
							     DDATA(v),
							     DDATA(penalty),
							     DDATA(fluxJacobian),
							     DDATA(fluxJacobian_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_LDG_upwind(PyObject* self,
                                                                PyObject* args)
{ 
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *u,
    *a,
    *phi_bc,
    *phi,
    *V,
    *penalty,
    *flux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &a,
                       &phi_bc,
                       &phi,
                       &V,
                       &penalty,
                       &flux))
    return NULL;
  if (ND(n) > 3)
    {
      calculateExteriorNumericalDiffusiveFlux_LDG_upwind(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(n)[1],
							 SHAPE(n)[2],
							 SHAPE(n)[3],
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 DDATA(n),
							 DDATA(u),
							 DDATA(a),
							 DDATA(phi_bc),
							 DDATA(phi),
							 DDATA(V),
							 DDATA(penalty),
							 DDATA(flux));
    }
  else
    {
      calculateGlobalExteriorNumericalDiffusiveFlux_LDG_upwind(SHAPE(exteriorElementBoundaries)[0],
							       SHAPE(V)[1],
							       SHAPE(V)[2],
							       SHAPE(V)[3],
							       IDATA(exteriorElementBoundaries),
							       IDATA(elementBoundaryElements),
							       IDATA(elementBoundaryLocalElementBoundaries),
							       DDATA(n),
							       DDATA(u),
							       DDATA(a),
							       DDATA(phi_bc),
							       DDATA(phi),
							       DDATA(V),
							       DDATA(penalty),
							       DDATA(flux));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_LDG_upwind_sd(PyObject* self,
								    PyObject* args)
{ 
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *u,
    *a,
    *phi_bc,
    *phi,
    *V,
    *penalty,
    *flux,
    *rowptr,
    *colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &a,
                       &phi_bc,
                       &phi,
                       &V,
                       &penalty,
                       &flux))
    return NULL;
  if (ND(n) > 3)
    {
      calculateExteriorNumericalDiffusiveFlux_LDG_upwind_sd(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(n)[1],
							    SHAPE(n)[2],
							    SHAPE(n)[3],
							    IDATA(rowptr),
							    IDATA(colind),
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    DDATA(n),
							    DDATA(u),
							    DDATA(a),
							    DDATA(phi_bc),
							    DDATA(phi),
							    DDATA(V),
							    DDATA(penalty),
							    DDATA(flux));
    }
  else
    {
      calculateGlobalExteriorNumericalDiffusiveFlux_LDG_upwind_sd(SHAPE(exteriorElementBoundaries)[0],
								  SHAPE(V)[1],
								  SHAPE(V)[2],
								  SHAPE(V)[3],
								  IDATA(rowptr),
								  IDATA(colind),
								  IDATA(exteriorElementBoundaries),
								  IDATA(elementBoundaryElements),
								  IDATA(elementBoundaryLocalElementBoundaries),
								  DDATA(n),
								  DDATA(u),
								  DDATA(a),
								  DDATA(phi_bc),
								  DDATA(phi),
								  DDATA(V),
								  DDATA(penalty),
								  DDATA(flux));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind(PyObject* self,
                                                                        PyObject* args)
{ 
  PyObject *isDiffusiveFluxBoundary,
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *a,
    *da,
    *dphi,
    *V,
    *dV,
    *dV_eb,
    *v,
    *penalty,
    *fluxJacobian,
    *fluxJacobian_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
		       &isDiffusiveFluxBoundary,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &dphi,
                       &V,
                       &dV,
                       &dV_eb,
                       &v,
                       &penalty,
                       &fluxJacobian,
                       &fluxJacobian_eb))
    return NULL;
  if (ND(n) > 3)
    {
      updateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind(IDATA(isDiffusiveFluxBoundary),
							      SHAPE(exteriorElementBoundaries)[0],
							      SHAPE(v)[1],
							      SHAPE(v)[2],
							      SHAPE(v)[3],
							      SHAPE(n)[3],
							      IDATA(exteriorElementBoundaries),
							      IDATA(elementBoundaryElements),
							      IDATA(elementBoundaryLocalElementBoundaries),
							      DDATA(n),
							      DDATA(a),
							      DDATA(da),
							      DDATA(dphi),
							      DDATA(V),
							      DDATA(dV),
							      DDATA(dV_eb),
							      DDATA(v),
							      DDATA(penalty),
							      DDATA(fluxJacobian),
							      DDATA(fluxJacobian_eb));
    }
  else
    {
      updateGlobalExteriorNumericalDiffusiveFluxJacobian_LDG_upwind(IDATA(isDiffusiveFluxBoundary),
								    SHAPE(exteriorElementBoundaries)[0],
								    SHAPE(dV)[1],
								    SHAPE(dV)[2],
								    SHAPE(dV)[3],
								    SHAPE(n)[2],
								    IDATA(exteriorElementBoundaries),
								    IDATA(elementBoundaryElements),
								    IDATA(elementBoundaryLocalElementBoundaries),
								    DDATA(n),
								    DDATA(a),
								    DDATA(da),
								    DDATA(dphi),
								    DDATA(V),
								    DDATA(dV),
								    DDATA(dV_eb),
								    DDATA(v),
								    DDATA(penalty),
								    DDATA(fluxJacobian),
								    DDATA(fluxJacobian_eb));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(PyObject* self,
									 PyObject* args)
{ 
  PyObject *isDiffusiveFluxBoundary,
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *a,
    *da,
    *dphi,
    *V,
    *dV,
    *dV_eb,
    *v,
    *penalty,
    *fluxJacobian,
    *fluxJacobian_eb,
    *rowptr,
    *colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &isDiffusiveFluxBoundary,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &a,
                       &da,
                       &dphi,
                       &V,
                       &dV,
                       &dV_eb,
                       &v,
                       &penalty,
                       &fluxJacobian,
                       &fluxJacobian_eb))
    return NULL;
  if (ND(n) > 3)
    {
      updateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(IDATA(isDiffusiveFluxBoundary),
								 SHAPE(exteriorElementBoundaries)[0],
								 SHAPE(v)[1],
								 SHAPE(v)[2],
								 SHAPE(v)[3],
								 SHAPE(n)[3],
								 IDATA(rowptr),
								 IDATA(colind),
								 IDATA(exteriorElementBoundaries),
								 IDATA(elementBoundaryElements),
								 IDATA(elementBoundaryLocalElementBoundaries),
								 DDATA(n),
								 DDATA(a),
								 DDATA(da),
								 DDATA(dphi),
								 DDATA(V),
								 DDATA(dV),
								 DDATA(dV_eb),
								 DDATA(v),
								 DDATA(penalty),
								 DDATA(fluxJacobian),
								 DDATA(fluxJacobian_eb));
    }
  else
    {
      updateGlobalExteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd(IDATA(isDiffusiveFluxBoundary),
								       SHAPE(exteriorElementBoundaries)[0],
								       SHAPE(dV)[1],
								       SHAPE(dV)[2],
								       SHAPE(dV)[3],
								       SHAPE(n)[2],
								       IDATA(rowptr),
								       IDATA(colind),
								       IDATA(exteriorElementBoundaries),
								       IDATA(elementBoundaryElements),
								       IDATA(elementBoundaryLocalElementBoundaries),
								       DDATA(n),
								       DDATA(a),
								       DDATA(da),
								       DDATA(dphi),
								       DDATA(V),
								       DDATA(dV),
								       DDATA(dV_eb),
								       DDATA(v),
								       DDATA(penalty),
								       DDATA(fluxJacobian),
								       DDATA(fluxJacobian_eb));
      
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateDiffusionMatrixSplittings_LDG_sd(PyObject* self,
							PyObject* args)
{ 
  PyObject *rowptr,
    *colind,
    *ebq_a,
    *q_a,
    *eb_aHat,
    *eb_aTilde,
    *aHat,
    *aTilde;
  int aSplit,nSpace;
  if(!PyArg_ParseTuple(args,"iiOOOOOOOO",
		       &aSplit,
		       &nSpace,
		       &rowptr,
		       &colind,
		       &ebq_a,
		       &q_a,
		       &eb_aHat,
		       &eb_aTilde,
		       &aHat,
		       &aTilde))

    return NULL;

  calculateDiffusionMatrixSplittings_LDG_sd(aSplit,
					    SHAPE(q_a)[0],
					    SHAPE(ebq_a)[1],
					    SHAPE(q_a)[1],
					    SHAPE(ebq_a)[2],
					    nSpace,
					    IDATA(rowptr),
					    IDATA(colind),
					    DDATA(ebq_a),
					    DDATA(q_a),
					    DDATA(eb_aHat),
					    DDATA(eb_aTilde),
					    DDATA(aHat),
					    DDATA(aTilde));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF(PyObject* self, 
							    PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_um,
    *n,
    *bc_f_m,  
    *bc_a_wm,      
    *bc_a_mw,      
    *bc_a_mm,      
    *bc_grad_phi_w,
    *bc_grad_phi_m,
    *bc_uw,        
    *bc_um,        
    *f_m,          
    *df_m_dw,      
    *a_wm,         
    *a_mw,         
    *a_mm,         
    *grad_phi_w,   
    *grad_phi_m,   
    *uw,           
    *um,           
    *penalty_w,    
    *penalty_m,
    *advectiveFlux_m,   
    *dadvectiveFlux_m_dw,   
    *diffusiveFlux_wm,
    *diffusiveFlux_mw,
    *diffusiveFlux_mm;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_um,
		       &n,
		       &bc_f_m,  
		       &bc_a_wm,      
		       &bc_a_mw,      
		       &bc_a_mm,      
		       &bc_grad_phi_w,
		       &bc_grad_phi_m,
		       &bc_uw,        
		       &bc_um,        
		       &f_m,          
		       &df_m_dw,      
		       &a_wm,         
		       &a_mw,         
		       &a_mm,         
		       &grad_phi_w,   
		       &grad_phi_m,   
		       &uw,           
		       &um,           
		       &penalty_w,    
		       &penalty_m,
		       &advectiveFlux_m,   
		       &dadvectiveFlux_m_dw,   
		       &diffusiveFlux_wm,
		       &diffusiveFlux_mw,
		       &diffusiveFlux_mm))
      return NULL;

    
    calculateGlobalExteriorNumericalFluxDarcyFCFF(SHAPE(exteriorElementBoundaries)[0],
						  SHAPE(n)[1],
						  SHAPE(n)[2],
						  IDATA(exteriorElementBoundaries),
						  IDATA(elementBoundaryElements),
						  IDATA(elementBoundaryLocalElementBoundaries),
						  IDATA(isDOFBoundary_uw),
						  IDATA(isDOFBoundary_um),
						  DDATA(n),
						  DDATA(bc_f_m),  
						  DDATA(bc_a_wm),      
						  DDATA(bc_a_mw),      
						  DDATA(bc_a_mm),      
						  DDATA(bc_grad_phi_w),
						  DDATA(bc_grad_phi_m),
						  DDATA(bc_uw),        
						  DDATA(bc_um),        
						  DDATA(f_m),          
						  DDATA(df_m_dw),      
						  DDATA(a_wm),         
						  DDATA(a_mw),         
						  DDATA(a_mm),         
						  DDATA(grad_phi_w),   
						  DDATA(grad_phi_m),   
						  DDATA(uw),           
						  DDATA(um),           
						  DDATA(penalty_w),    
						  DDATA(penalty_m),
						  DDATA(advectiveFlux_m),   
						  DDATA(dadvectiveFlux_m_dw),   
						  DDATA(diffusiveFlux_wm),
						  DDATA(diffusiveFlux_mw),
						  DDATA(diffusiveFlux_mm));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF_sd(PyObject* self, 
							       PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_um,
    *n,
    *bc_f_m,  
    *bc_a_wm,      
    *bc_a_mw,      
    *bc_a_mm,      
    *bc_grad_phi_w,
    *bc_grad_phi_m,
    *bc_uw,        
    *bc_um,        
    *f_m,          
    *df_m_dw,      
    *a_wm,         
    *a_mw,         
    *a_mm,         
    *grad_phi_w,   
    *grad_phi_m,   
    *uw,           
    *um,           
    *penalty_w,    
    *penalty_m,
    *advectiveFlux_m,   
    *dadvectiveFlux_m_dw,   
    *diffusiveFlux_wm,
    *diffusiveFlux_mw,
    *diffusiveFlux_mm,
    *rowptr_wm,
    *colind_wm,
    *rowptr_mw,
    *colind_mw,
    *rowptr_mm,
    *colind_mm;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &rowptr_wm,
		       &colind_wm,
		       &rowptr_mw,
		       &colind_mw,
		       &rowptr_mm,
		       &colind_mm,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_um,
		       &n,
		       &bc_f_m,  
		       &bc_a_wm,      
		       &bc_a_mw,      
		       &bc_a_mm,      
		       &bc_grad_phi_w,
		       &bc_grad_phi_m,
		       &bc_uw,        
		       &bc_um,        
		       &f_m,          
		       &df_m_dw,      
		       &a_wm,         
		       &a_mw,         
		       &a_mm,         
		       &grad_phi_w,   
		       &grad_phi_m,   
		       &uw,           
		       &um,           
		       &penalty_w,    
		       &penalty_m,
		       &advectiveFlux_m,   
		       &dadvectiveFlux_m_dw,   
		       &diffusiveFlux_wm,
		       &diffusiveFlux_mw,
		       &diffusiveFlux_mm))
      return NULL;

    
  calculateGlobalExteriorNumericalFluxDarcyFCFF_sd(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(n)[1],
						   SHAPE(n)[2],
						   IDATA(rowptr_wm),
						   IDATA(colind_wm),
						   IDATA(rowptr_mw),
						   IDATA(colind_mw),
						   IDATA(rowptr_mm),
						   IDATA(colind_mm),
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary_uw),
						   IDATA(isDOFBoundary_um),
						   DDATA(n),
						   DDATA(bc_f_m),  
						   DDATA(bc_a_wm),      
						   DDATA(bc_a_mw),      
						   DDATA(bc_a_mm),      
						   DDATA(bc_grad_phi_w),
						   DDATA(bc_grad_phi_m),
						   DDATA(bc_uw),        
						   DDATA(bc_um),        
						   DDATA(f_m),          
						   DDATA(df_m_dw),      
						   DDATA(a_wm),         
						   DDATA(a_mw),         
						   DDATA(a_mm),         
						   DDATA(grad_phi_w),   
						   DDATA(grad_phi_m),   
						   DDATA(uw),           
						   DDATA(um),           
						   DDATA(penalty_w),    
						   DDATA(penalty_m),
						   DDATA(advectiveFlux_m),   
						   DDATA(dadvectiveFlux_m_dw),   
						   DDATA(diffusiveFlux_wm),
						   DDATA(diffusiveFlux_mw),
						   DDATA(diffusiveFlux_mm));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian(PyObject* self, 
										  PyObject* args)
{
  PyObject *l2g, /*for now assumes both solution spaces are the same!*/
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_um,
    *n,
    *f_m, 
    *df_m_dw,
    *a_wm,   
    *da_wm_dw,
    *da_wm_dm,
    *a_mw,    
    *da_mw_dw,
    *da_mw_dm,
    *a_mm,    
    *da_mm_dw,
    *da_mm_dm,
    *grad_phi_w,
    *grad_phi_m,
    *dphi_w_w,  
    *dphi_w_m,  
    *dphi_m_w,  
    *dphi_m_m,  
    *u_w,        
    *u_m,        
    *v,         
    *grad_v,    
    *penalty_w,    
    *penalty_m,
    *fluxJacobian_ww,
    *fluxJacobian_wm,
    *fluxJacobian_mw,
    *fluxJacobian_mm;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &l2g, 
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_um,
		       &n,
		       &f_m, 
		       &df_m_dw,
		       &a_wm,   
		       &da_wm_dw,
		       &da_wm_dm,
		       &a_mw,    
		       &da_mw_dw,
		       &da_mw_dm,
		       &a_mm,    
		       &da_mm_dw,
		       &da_mm_dm,
		       &grad_phi_w,
		       &grad_phi_m,
		       &dphi_w_w,  
		       &dphi_w_m,  
		       &dphi_m_w,  
		       &dphi_m_m,  
		       &u_w,        
		       &u_m,        
		       &v,         
		       &grad_v,    
		       &penalty_w,    
		       &penalty_m,
		       &fluxJacobian_ww,
		       &fluxJacobian_wm,
		       &fluxJacobian_mw,
		       &fluxJacobian_mm))
    return NULL;

  calculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
								      SHAPE(n)[1],
								      SHAPE(n)[2],
								      SHAPE(l2g)[1],
								      IDATA(l2g), 
								      IDATA(exteriorElementBoundaries),
								      IDATA(elementBoundaryElements),
								      IDATA(elementBoundaryLocalElementBoundaries),
								      IDATA(isDOFBoundary_uw),
								      IDATA(isDOFBoundary_um),
								      DDATA(n),
								      DDATA(f_m), 
								      DDATA(df_m_dw),
								      DDATA(a_wm),   
								      DDATA(da_wm_dw),
								      DDATA(da_wm_dm),
								      DDATA(a_mw),    
								      DDATA(da_mw_dw),
								      DDATA(da_mw_dm),
								      DDATA(a_mm),    
								      DDATA(da_mm_dw),
								      DDATA(da_mm_dm),
								      DDATA(grad_phi_w),
								      DDATA(grad_phi_m),
								      DDATA(dphi_w_w),  
								      DDATA(dphi_w_m),  
								      DDATA(dphi_m_w),  
								      DDATA(dphi_m_m),  
								      DDATA(u_w),        
								      DDATA(u_m),        
								      DDATA(v),         
								      DDATA(grad_v),    
								      DDATA(penalty_w),    
								      DDATA(penalty_m),
								      DDATA(fluxJacobian_ww),
								      DDATA(fluxJacobian_wm),
								      DDATA(fluxJacobian_mw),
								      DDATA(fluxJacobian_mm));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian_sd(PyObject* self, 
										     PyObject* args)
{
  PyObject *l2g, /*for now assumes both solution spaces are the same!*/
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_um,
    *n,
    *f_m, 
    *df_m_dw,
    *a_wm,   
    *da_wm_dw,
    *da_wm_dm,
    *a_mw,    
    *da_mw_dw,
    *da_mw_dm,
    *a_mm,    
    *da_mm_dw,
    *da_mm_dm,
    *grad_phi_w,
    *grad_phi_m,
    *dphi_w_w,  
    *dphi_w_m,  
    *dphi_m_w,  
    *dphi_m_m,  
    *u_w,        
    *u_m,        
    *v,         
    *grad_v,    
    *penalty_w,    
    *penalty_m,
    *fluxJacobian_ww,
    *fluxJacobian_wm,
    *fluxJacobian_mw,
    *fluxJacobian_mm,
    *rowptr_wm,
    *colind_wm,
    *rowptr_mw,
    *colind_mw,
    *rowptr_mm,
    *colind_mm;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &rowptr_wm,
		       &colind_wm,
		       &rowptr_mw,
		       &colind_mw,
		       &rowptr_mm,
		       &colind_mm,
		       &l2g, 
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_um,
		       &n,
		       &f_m, 
		       &df_m_dw,
		       &a_wm,   
		       &da_wm_dw,
		       &da_wm_dm,
		       &a_mw,    
		       &da_mw_dw,
		       &da_mw_dm,
		       &a_mm,    
		       &da_mm_dw,
		       &da_mm_dm,
		       &grad_phi_w,
		       &grad_phi_m,
		       &dphi_w_w,  
		       &dphi_w_m,  
		       &dphi_m_w,  
		       &dphi_m_m,  
		       &u_w,        
		       &u_m,        
		       &v,         
		       &grad_v,    
		       &penalty_w,    
		       &penalty_m,
		       &fluxJacobian_ww,
		       &fluxJacobian_wm,
		       &fluxJacobian_mw,
		       &fluxJacobian_mm))
    return NULL;

  calculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian_sd(SHAPE(exteriorElementBoundaries)[0],
								      SHAPE(n)[1],
								      SHAPE(n)[2],
								      SHAPE(l2g)[1],
									 IDATA(rowptr_wm),
									 IDATA(colind_wm),
									 IDATA(rowptr_mw),
									 IDATA(colind_mw),
									 IDATA(rowptr_mm),
									 IDATA(colind_mm),
								      IDATA(l2g), 
								      IDATA(exteriorElementBoundaries),
								      IDATA(elementBoundaryElements),
								      IDATA(elementBoundaryLocalElementBoundaries),
								      IDATA(isDOFBoundary_uw),
								      IDATA(isDOFBoundary_um),
								      DDATA(n),
								      DDATA(f_m), 
								      DDATA(df_m_dw),
								      DDATA(a_wm),   
								      DDATA(da_wm_dw),
								      DDATA(da_wm_dm),
								      DDATA(a_mw),    
								      DDATA(da_mw_dw),
								      DDATA(da_mw_dm),
								      DDATA(a_mm),    
								      DDATA(da_mm_dw),
								      DDATA(da_mm_dm),
								      DDATA(grad_phi_w),
								      DDATA(grad_phi_m),
								      DDATA(dphi_w_w),  
								      DDATA(dphi_w_m),  
								      DDATA(dphi_m_w),  
								      DDATA(dphi_m_m),  
								      DDATA(u_w),        
								      DDATA(u_m),        
								      DDATA(v),         
								      DDATA(grad_v),    
								      DDATA(penalty_w),    
								      DDATA(penalty_m),
								      DDATA(fluxJacobian_ww),
								      DDATA(fluxJacobian_wm),
								      DDATA(fluxJacobian_mw),
								      DDATA(fluxJacobian_mm));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC(PyObject* self, 
							  PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *bc_a_ww,      
    *bc_a_nn,      
    *bc_grad_phi_w,
    *bc_grad_phi_n,
    *bc_uw,        
    *bc_un,
    *bc_psi_n,
    *a_ww,         
    *a_nn,         
    *grad_phi_w,   
    *grad_phi_n,   
    *uw,           
    *un,           
    *psi_n,
    *penalty_w,    
    *penalty_n,
    *diffusiveFlux_ww,
    *diffusiveFlux_nn;

  if(!PyArg_ParseTuple(args,"OOOOOiiOOOOOOOOOOOOOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &bc_a_ww,      
		       &bc_a_nn,      
		       &bc_grad_phi_w,
		       &bc_grad_phi_n,
		       &bc_uw,        
		       &bc_un,   
		       &bc_psi_n,
		       &a_ww,         
		       &a_nn,         
		       &grad_phi_w,   
		       &grad_phi_n,   
		       &uw,           
		       &un,      
		       &psi_n,
		       &penalty_w,    
		       &penalty_n,
		       &diffusiveFlux_ww,
		       &diffusiveFlux_nn))
      return NULL;

    
    calculateGlobalExteriorNumericalFluxDarcyFC(SHAPE(exteriorElementBoundaries)[0],
						SHAPE(n)[1],
						SHAPE(n)[2],
						IDATA(exteriorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						IDATA(isDOFBoundary_uw),
						IDATA(isDOFBoundary_un),
						fluxBoundaryFlags_uw,
						fluxBoundaryFlags_un,
						DDATA(n),
						DDATA(bc_a_ww),      
						DDATA(bc_a_nn),      
						DDATA(bc_grad_phi_w),
						DDATA(bc_grad_phi_n),
						DDATA(bc_uw),        
						DDATA(bc_un),        
						DDATA(bc_psi_n),
						DDATA(a_ww),         
						DDATA(a_nn),         
						DDATA(grad_phi_w),   
						DDATA(grad_phi_n),   
						DDATA(uw),           
						DDATA(un),    
						DDATA(psi_n),
						DDATA(penalty_w),    
						DDATA(penalty_n),
						DDATA(diffusiveFlux_ww),
						DDATA(diffusiveFlux_nn));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC_sd(PyObject* self, 
							     PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *bc_a_ww,      
    *bc_a_nn,      
    *bc_grad_phi_w,
    *bc_grad_phi_n,
    *bc_uw,        
    *bc_un,
    *bc_psi_n,
    *a_ww,         
    *a_nn,         
    *grad_phi_w,   
    *grad_phi_n,   
    *uw,           
    *un,           
    *psi_n,
    *penalty_w,    
    *penalty_n,
    *diffusiveFlux_ww,
    *diffusiveFlux_nn,
    *rowptr_ww,
    *colind_ww,
    *rowptr_nn,
    *colind_nn;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOiiOOOOOOOOOOOOOOOOOOO",
		       &rowptr_ww,
		       &colind_ww,
		       &rowptr_nn,
		       &colind_nn,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &bc_a_ww,      
		       &bc_a_nn,      
		       &bc_grad_phi_w,
		       &bc_grad_phi_n,
		       &bc_uw,        
		       &bc_un,   
		       &bc_psi_n,
		       &a_ww,         
		       &a_nn,         
		       &grad_phi_w,   
		       &grad_phi_n,   
		       &uw,           
		       &un,      
		       &psi_n,
		       &penalty_w,    
		       &penalty_n,
		       &diffusiveFlux_ww,
		       &diffusiveFlux_nn))
      return NULL;

    
  calculateGlobalExteriorNumericalFluxDarcyFC_sd(SHAPE(exteriorElementBoundaries)[0],
						 SHAPE(n)[1],
						 SHAPE(n)[2],
						 IDATA(rowptr_ww),
						 IDATA(colind_ww),
						 IDATA(rowptr_nn),
						 IDATA(colind_nn),
						 IDATA(exteriorElementBoundaries),
						 IDATA(elementBoundaryElements),
						 IDATA(elementBoundaryLocalElementBoundaries),
						 IDATA(isDOFBoundary_uw),
						 IDATA(isDOFBoundary_un),
						 fluxBoundaryFlags_uw,
						 fluxBoundaryFlags_un,
						 DDATA(n),
						 DDATA(bc_a_ww),      
						 DDATA(bc_a_nn),      
						 DDATA(bc_grad_phi_w),
						 DDATA(bc_grad_phi_n),
						 DDATA(bc_uw),        
						 DDATA(bc_un),        
						 DDATA(bc_psi_n),
						 DDATA(a_ww),         
						 DDATA(a_nn),         
						 DDATA(grad_phi_w),   
						 DDATA(grad_phi_n),   
						 DDATA(uw),           
						 DDATA(un),    
						 DDATA(psi_n),
						 DDATA(penalty_w),    
						 DDATA(penalty_n),
						 DDATA(diffusiveFlux_ww),
						 DDATA(diffusiveFlux_nn));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian(PyObject* self, 
										PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *l2g, /*for now assumes both solution spaces are the same!*/
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *a_ww,   
    *da_ww_dw,
    *da_ww_dn,
    *a_nn,    
    *da_nn_dw,
    *da_nn_dn,
    *grad_phi_w,
    *grad_phi_n,
    *dphi_w_w,  
    *dphi_w_n,  
    *dphi_n_w,  
    *dphi_n_n,  
    *u_w,        
    *u_n,        
    *psi_n,
    *dpsi_n_dsw,
    *dpsi_n_dpsiw,
    *v,         
    *grad_v,    
    *penalty_w,    
    *penalty_n,
    *fluxJacobian_ww,
    *fluxJacobian_wn,
    *fluxJacobian_nw,
    *fluxJacobian_nn;

  if(!PyArg_ParseTuple(args,"OOOOOOiiOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &l2g, 
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &a_ww,   
		       &da_ww_dw,
		       &da_ww_dn,
		       &a_nn,    
		       &da_nn_dw,
		       &da_nn_dn,
		       &grad_phi_w,
		       &grad_phi_n,
		       &dphi_w_w,  
		       &dphi_w_n,  
		       &dphi_n_w,  
		       &dphi_n_n,  
		       &u_w,        
		       &u_n,        
		       &psi_n,
		       &dpsi_n_dsw,
		       &dpsi_n_dpsiw,
		       &v,         
		       &grad_v,    
		       &penalty_w,    
		       &penalty_n,
		       &fluxJacobian_ww,
		       &fluxJacobian_wn,
		       &fluxJacobian_nw,
		       &fluxJacobian_nn))
    return NULL;

  calculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
								    SHAPE(n)[1],
								    SHAPE(n)[2],
								    SHAPE(l2g)[1],
								    IDATA(l2g), 
								    IDATA(exteriorElementBoundaries),
								    IDATA(elementBoundaryElements),
								    IDATA(elementBoundaryLocalElementBoundaries),
								    IDATA(isDOFBoundary_uw),
								    IDATA(isDOFBoundary_un),
								    fluxBoundaryFlags_uw,
								    fluxBoundaryFlags_un,
								    DDATA(n),
								    DDATA(a_ww),   
								    DDATA(da_ww_dw),
								    DDATA(da_ww_dn),
								    DDATA(a_nn),    
								    DDATA(da_nn_dw),
								    DDATA(da_nn_dn),
								    DDATA(grad_phi_w),
								    DDATA(grad_phi_n),
								    DDATA(dphi_w_w),  
								    DDATA(dphi_w_n),  
								    DDATA(dphi_n_w),  
								    DDATA(dphi_n_n),  
								    DDATA(u_w),        
								    DDATA(u_n),       
								    DDATA(psi_n),
								    DDATA(dpsi_n_dsw),
								    DDATA(dpsi_n_dpsiw),
								    DDATA(v),         
								    DDATA(grad_v),    
								    DDATA(penalty_w),    
								    DDATA(penalty_n),
								    DDATA(fluxJacobian_ww),
								    DDATA(fluxJacobian_wn),
								    DDATA(fluxJacobian_nw),
								    DDATA(fluxJacobian_nn));

  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian_sd(PyObject* self, 
										   PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *l2g, /*for now assumes both solution spaces are the same!*/
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *a_ww,   
    *da_ww_dw,
    *da_ww_dn,
    *a_nn,    
    *da_nn_dw,
    *da_nn_dn,
    *grad_phi_w,
    *grad_phi_n,
    *dphi_w_w,  
    *dphi_w_n,  
    *dphi_n_w,  
    *dphi_n_n,  
    *u_w,        
    *u_n,        
    *psi_n,
    *dpsi_n_dsw,
    *dpsi_n_dpsiw,
    *v,         
    *grad_v,    
    *penalty_w,    
    *penalty_n,
    *fluxJacobian_ww,
    *fluxJacobian_wn,
    *fluxJacobian_nw,
    *fluxJacobian_nn,
    *rowptr_ww,
    *colind_ww,
    *rowptr_nn,
    *colind_nn;
  
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOiiOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &rowptr_ww,
		       &colind_ww,
		       &rowptr_nn,
		       &colind_nn,
		       &l2g, 
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &a_ww,   
		       &da_ww_dw,
		       &da_ww_dn,
		       &a_nn,    
		       &da_nn_dw,
		       &da_nn_dn,
		       &grad_phi_w,
		       &grad_phi_n,
		       &dphi_w_w,  
		       &dphi_w_n,  
		       &dphi_n_w,  
		       &dphi_n_n,  
		       &u_w,        
		       &u_n,        
		       &psi_n,
		       &dpsi_n_dsw,
		       &dpsi_n_dpsiw,
		       &v,         
		       &grad_v,    
		       &penalty_w,    
		       &penalty_n,
		       &fluxJacobian_ww,
		       &fluxJacobian_wn,
		       &fluxJacobian_nw,
		       &fluxJacobian_nn))
    return NULL;

  calculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian_sd(SHAPE(exteriorElementBoundaries)[0],
								       SHAPE(n)[1],
								       SHAPE(n)[2],
								       SHAPE(l2g)[1],
								       IDATA(rowptr_ww),
								       IDATA(colind_ww),
								       IDATA(rowptr_nn),
								       IDATA(colind_nn),
								       IDATA(l2g), 
								       IDATA(exteriorElementBoundaries),
								       IDATA(elementBoundaryElements),
								       IDATA(elementBoundaryLocalElementBoundaries),
								       IDATA(isDOFBoundary_uw),
								       IDATA(isDOFBoundary_un),
								       fluxBoundaryFlags_uw,
								       fluxBoundaryFlags_un,
								       DDATA(n),
								       DDATA(a_ww),   
								       DDATA(da_ww_dw),
								       DDATA(da_ww_dn),
								       DDATA(a_nn),    
								       DDATA(da_nn_dw),
								       DDATA(da_nn_dn),
								       DDATA(grad_phi_w),
								       DDATA(grad_phi_n),
								       DDATA(dphi_w_w),  
								       DDATA(dphi_w_n),  
								       DDATA(dphi_n_w),  
								       DDATA(dphi_n_n),  
								       DDATA(u_w),        
								       DDATA(u_n),       
								       DDATA(psi_n),
								       DDATA(dpsi_n_dsw),
								       DDATA(dpsi_n_dpsiw),
								       DDATA(v),         
								       DDATA(grad_v),    
								       DDATA(penalty_w),    
								       DDATA(penalty_n),
								       DDATA(fluxJacobian_ww),
								       DDATA(fluxJacobian_wn),
								       DDATA(fluxJacobian_nw),
								       DDATA(fluxJacobian_nn));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalAdvectiveFlux_DarcyFC(PyObject* self, 
								    PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*n,*bc_sw,*bc_psiw,
    *bc_fw,*bc_dfw_dsw,*bc_dfw_dpsiw,*bc_fn,*bc_dfn_dsw,*bc_dfn_dpsiw,*sw,*psiw,*fw,*dfw_dsw,*dfw_dpsiw,*fn,*dfn_dsw,*dfn_dpsiw,
    *advectiveFlux_w,*dadvectiveFlux_w_dsw_left,*dadvectiveFlux_w_dpsiw_left,
    *advectiveFlux_n,*dadvectiveFlux_n_dsw_left,*dadvectiveFlux_n_dpsiw_left,
    *isDOFBoundary_sw,*isDOFBoundary_psiw;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary_sw,
		       &isDOFBoundary_psiw,
                       &n,
                       &bc_sw,
		       &bc_psiw,
                       &bc_fw,
		       &bc_dfw_dsw,
                       &bc_dfw_dpsiw,
		       &bc_fn,
		       &bc_dfn_dsw,
                       &bc_dfn_dpsiw,
                       &sw,
		       &psiw,
                       &fw,
                       &dfw_dsw,
		       &dfw_dpsiw,
		       &fn,
		       &dfn_dsw,
		       &dfn_dpsiw,
                       &advectiveFlux_w,
                       &dadvectiveFlux_w_dsw_left,
                       &dadvectiveFlux_w_dpsiw_left,
                       &advectiveFlux_n,
                       &dadvectiveFlux_n_dsw_left,
                       &dadvectiveFlux_n_dpsiw_left))
    return NULL;
  assert(ND(fw) == 3);
  calculateGlobalExteriorNumericalAdvectiveFlux_DarcyFC(SHAPE(exteriorElementBoundaries)[0],
							SHAPE(fw)[1],
							SHAPE(fw)[2],
							IDATA(exteriorElementBoundaries),
							IDATA(elementBoundaryElements),
							IDATA(elementBoundaryLocalElementBoundaries),
							IDATA(isDOFBoundary_sw),
							IDATA(isDOFBoundary_psiw),
							DDATA(n),
							DDATA(bc_sw),
							DDATA(bc_psiw),
							DDATA(bc_fw),
							DDATA(bc_dfw_dsw),
							DDATA(bc_dfw_dpsiw),
							DDATA(bc_fn),
							DDATA(bc_dfn_dsw),
							DDATA(bc_dfn_dpsiw),
							DDATA(sw),
							DDATA(psiw),
							DDATA(fw),
							DDATA(dfw_dsw),
							DDATA(dfw_dpsiw),
							DDATA(fn),
							DDATA(dfn_dsw),
							DDATA(dfn_dpsiw),
							DDATA(advectiveFlux_w),
							DDATA(dadvectiveFlux_w_dsw_left),
							DDATA(dadvectiveFlux_w_dpsiw_left),
							DDATA(advectiveFlux_n),
							DDATA(dadvectiveFlux_n_dsw_left),
							DDATA(dadvectiveFlux_n_dpsiw_left));
    
  Py_INCREF(Py_None);
  return Py_None;
}


/*start FCPP routines */

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP(PyObject* self, 
							    PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *bc_a_ww,      
    *bc_a_nn,      
    *bc_grad_phi_w,
    *bc_grad_phi_n,
    *bc_uw,        
    *bc_un,
    *bc_psi_n,
    *a_ww,         
    *a_nn,         
    *grad_phi_w,   
    *grad_phi_n,   
    *uw,           
    *un,           
    *psi_n,
    *penalty_w,    
    *penalty_n,
    *diffusiveFlux_ww,
    *diffusiveFlux_nn;

  if(!PyArg_ParseTuple(args,"OOOOOiiOOOOOOOOOOOOOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &bc_a_ww,      
		       &bc_a_nn,      
		       &bc_grad_phi_w,
		       &bc_grad_phi_n,
		       &bc_uw,        
		       &bc_un,   
		       &bc_psi_n,
		       &a_ww,         
		       &a_nn,         
		       &grad_phi_w,   
		       &grad_phi_n,   
		       &uw,           
		       &un,      
		       &psi_n,
		       &penalty_w,    
		       &penalty_n,
		       &diffusiveFlux_ww,
		       &diffusiveFlux_nn))
      return NULL;

    
    calculateGlobalExteriorNumericalFluxDarcyFCPP(SHAPE(exteriorElementBoundaries)[0],
						  SHAPE(n)[1],
						  SHAPE(n)[2],
						  IDATA(exteriorElementBoundaries),
						  IDATA(elementBoundaryElements),
						  IDATA(elementBoundaryLocalElementBoundaries),
						  IDATA(isDOFBoundary_uw),
						  IDATA(isDOFBoundary_un),
						  fluxBoundaryFlags_uw,
						  fluxBoundaryFlags_un,
						  DDATA(n),
						  DDATA(bc_a_ww),      
						  DDATA(bc_a_nn),      
						  DDATA(bc_grad_phi_w),
						  DDATA(bc_grad_phi_n),
						  DDATA(bc_uw),        
						  DDATA(bc_un),        
						  DDATA(bc_psi_n),
						  DDATA(a_ww),         
						  DDATA(a_nn),         
						  DDATA(grad_phi_w),   
						  DDATA(grad_phi_n),   
						  DDATA(uw),           
						  DDATA(un),    
						  DDATA(psi_n),
						  DDATA(penalty_w),    
						  DDATA(penalty_n),
						  DDATA(diffusiveFlux_ww),
						  DDATA(diffusiveFlux_nn));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP_sd(PyObject* self, 
							       PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *bc_a_ww,      
    *bc_a_nn,      
    *bc_grad_phi_w,
    *bc_grad_phi_n,
    *bc_uw,        
    *bc_un,
    *bc_psi_n,
    *a_ww,         
    *a_nn,         
    *grad_phi_w,   
    *grad_phi_n,   
    *uw,           
    *un,           
    *psi_n,
    *penalty_w,    
    *penalty_n,
    *diffusiveFlux_ww,
    *diffusiveFlux_nn,
    *rowptr_ww,
    *colind_ww,
    *rowptr_nn,
    *colind_nn;

  if(!PyArg_ParseTuple(args,"OOOOOOOOOiiOOOOOOOOOOOOOOOOOOO",
		       &rowptr_ww,
		       &colind_ww,
		       &rowptr_nn,
		       &colind_nn,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &bc_a_ww,      
		       &bc_a_nn,      
		       &bc_grad_phi_w,
		       &bc_grad_phi_n,
		       &bc_uw,        
		       &bc_un,   
		       &bc_psi_n,
		       &a_ww,         
		       &a_nn,         
		       &grad_phi_w,   
		       &grad_phi_n,   
		       &uw,           
		       &un,      
		       &psi_n,
		       &penalty_w,    
		       &penalty_n,
		       &diffusiveFlux_ww,
		       &diffusiveFlux_nn))
      return NULL;

    
  calculateGlobalExteriorNumericalFluxDarcyFCPP_sd(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(n)[1],
						   SHAPE(n)[2],
						   IDATA(rowptr_ww),
						   IDATA(colind_ww),
						   IDATA(rowptr_nn),
						   IDATA(colind_nn),
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   IDATA(isDOFBoundary_uw),
						   IDATA(isDOFBoundary_un),
						   fluxBoundaryFlags_uw,
						   fluxBoundaryFlags_un,
						   DDATA(n),
						   DDATA(bc_a_ww),      
						   DDATA(bc_a_nn),      
						   DDATA(bc_grad_phi_w),
						   DDATA(bc_grad_phi_n),
						   DDATA(bc_uw),        
						   DDATA(bc_un),        
						   DDATA(bc_psi_n),
						   DDATA(a_ww),         
						   DDATA(a_nn),         
						   DDATA(grad_phi_w),   
						   DDATA(grad_phi_n),   
						   DDATA(uw),           
						   DDATA(un),    
						   DDATA(psi_n),
						   DDATA(penalty_w),    
						   DDATA(penalty_n),
						   DDATA(diffusiveFlux_ww),
						   DDATA(diffusiveFlux_nn));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian(PyObject* self, 
										  PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *l2g, /*for now assumes both solution spaces are the same!*/
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *a_ww,   
    *da_ww_dw,
    *da_ww_dn,
    *a_nn,    
    *da_nn_dw,
    *da_nn_dn,
    *grad_phi_w,
    *grad_phi_n,
    *dphi_w_w,  
    *dphi_w_n,  
    *dphi_n_w,  
    *dphi_n_n,  
    *u_w,        
    *u_n,        
    *psi_n,
    *dpsi_n_dpsiw,
    *dpsi_n_dpsic,
    *v,         
    *grad_v,    
    *penalty_w,    
    *penalty_n,
    *fluxJacobian_ww,
    *fluxJacobian_wn,
    *fluxJacobian_nw,
    *fluxJacobian_nn;

  if(!PyArg_ParseTuple(args,"OOOOOOiiOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &l2g, 
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &a_ww,   
		       &da_ww_dw,
		       &da_ww_dn,
		       &a_nn,    
		       &da_nn_dw,
		       &da_nn_dn,
		       &grad_phi_w,
		       &grad_phi_n,
		       &dphi_w_w,  
		       &dphi_w_n,  
		       &dphi_n_w,  
		       &dphi_n_n,  
		       &u_w,        
		       &u_n,        
		       &psi_n,
		       &dpsi_n_dpsiw,
		       &dpsi_n_dpsic,
		       &v,         
		       &grad_v,    
		       &penalty_w,    
		       &penalty_n,
		       &fluxJacobian_ww,
		       &fluxJacobian_wn,
		       &fluxJacobian_nw,
		       &fluxJacobian_nn))
    return NULL;

  calculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian(SHAPE(exteriorElementBoundaries)[0],
								      SHAPE(n)[1],
								      SHAPE(n)[2],
								      SHAPE(l2g)[1],
								      IDATA(l2g), 
								      IDATA(exteriorElementBoundaries),
								      IDATA(elementBoundaryElements),
								      IDATA(elementBoundaryLocalElementBoundaries),
								      IDATA(isDOFBoundary_uw),
								      IDATA(isDOFBoundary_un),
								      fluxBoundaryFlags_uw,
								      fluxBoundaryFlags_un,
								      DDATA(n),
								      DDATA(a_ww),   
								      DDATA(da_ww_dw),
								      DDATA(da_ww_dn),
								      DDATA(a_nn),    
								      DDATA(da_nn_dw),
								      DDATA(da_nn_dn),
								      DDATA(grad_phi_w),
								      DDATA(grad_phi_n),
								      DDATA(dphi_w_w),  
								      DDATA(dphi_w_n),  
								      DDATA(dphi_n_w),  
								      DDATA(dphi_n_n),  
								      DDATA(u_w),        
								      DDATA(u_n),       
								      DDATA(psi_n),
								      DDATA(dpsi_n_dpsiw),
								      DDATA(dpsi_n_dpsic),
								      DDATA(v),         
								      DDATA(grad_v),    
								      DDATA(penalty_w),    
								      DDATA(penalty_n),
								      DDATA(fluxJacobian_ww),
								      DDATA(fluxJacobian_wn),
								      DDATA(fluxJacobian_nw),
								      DDATA(fluxJacobian_nn));
  
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian_sd(PyObject* self, 
										     PyObject* args)
{
  int fluxBoundaryFlags_uw,fluxBoundaryFlags_un;
  PyObject *l2g, /*for now assumes both solution spaces are the same!*/
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_uw,
    *isDOFBoundary_un,
    *n,
    *a_ww,   
    *da_ww_dw,
    *da_ww_dn,
    *a_nn,    
    *da_nn_dw,
    *da_nn_dn,
    *grad_phi_w,
    *grad_phi_n,
    *dphi_w_w,  
    *dphi_w_n,  
    *dphi_n_w,  
    *dphi_n_n,  
    *u_w,        
    *u_n,        
    *psi_n,
    *dpsi_n_dpsiw,
    *dpsi_n_dpsic,
    *v,         
    *grad_v,    
    *penalty_w,    
    *penalty_n,
    *fluxJacobian_ww,
    *fluxJacobian_wn,
    *fluxJacobian_nw,
    *fluxJacobian_nn,
    *rowptr_ww,
    *colind_ww,
    *rowptr_nn,
    *colind_nn;
  
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOiiOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &rowptr_ww,
		       &colind_ww,
		       &rowptr_nn,
		       &colind_nn,
		       &l2g, 
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_uw,
		       &isDOFBoundary_un,
		       &fluxBoundaryFlags_uw,
		       &fluxBoundaryFlags_un,
		       &n,
		       &a_ww,   
		       &da_ww_dw,
		       &da_ww_dn,
		       &a_nn,    
		       &da_nn_dw,
		       &da_nn_dn,
		       &grad_phi_w,
		       &grad_phi_n,
		       &dphi_w_w,  
		       &dphi_w_n,  
		       &dphi_n_w,  
		       &dphi_n_n,  
		       &u_w,        
		       &u_n,        
		       &psi_n,
		       &dpsi_n_dpsiw,
		       &dpsi_n_dpsic,
		       &v,         
		       &grad_v,    
		       &penalty_w,    
		       &penalty_n,
		       &fluxJacobian_ww,
		       &fluxJacobian_wn,
		       &fluxJacobian_nw,
		       &fluxJacobian_nn))
    return NULL;

  calculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian_sd(SHAPE(exteriorElementBoundaries)[0],
									 SHAPE(n)[1],
									 SHAPE(n)[2],
									 SHAPE(l2g)[1],
									 IDATA(rowptr_ww),
									 IDATA(colind_ww),
									 IDATA(rowptr_nn),
									 IDATA(colind_nn),
									 IDATA(l2g), 
									 IDATA(exteriorElementBoundaries),
									 IDATA(elementBoundaryElements),
									 IDATA(elementBoundaryLocalElementBoundaries),
									 IDATA(isDOFBoundary_uw),
									 IDATA(isDOFBoundary_un),
									 fluxBoundaryFlags_uw,
									 fluxBoundaryFlags_un,
									 DDATA(n),
									 DDATA(a_ww),   
									 DDATA(da_ww_dw),
									 DDATA(da_ww_dn),
									 DDATA(a_nn),    
									 DDATA(da_nn_dw),
									 DDATA(da_nn_dn),
									 DDATA(grad_phi_w),
									 DDATA(grad_phi_n),
									 DDATA(dphi_w_w),  
									 DDATA(dphi_w_n),  
									 DDATA(dphi_n_w),  
									 DDATA(dphi_n_n),  
									 DDATA(u_w),        
									 DDATA(u_n),       
									 DDATA(psi_n),
									 DDATA(dpsi_n_dpsiw),
									 DDATA(dpsi_n_dpsic),
									 DDATA(v),         
									 DDATA(grad_v),    
									 DDATA(penalty_w),    
									 DDATA(penalty_n),
									 DDATA(fluxJacobian_ww),
									 DDATA(fluxJacobian_wn),
									 DDATA(fluxJacobian_nw),
									 DDATA(fluxJacobian_nn));
  
  Py_INCREF(Py_None);
  return Py_None;
}


/*end FCPP routines*/
static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcySplitPressure(PyObject* self, 
								     PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_psiw,*bc_psin,
    *a,*grad_phi,*psiw,*psin,*penalty,*diffusiveFlux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_psiw,
		       &bc_psin,
                       &a,
                       &grad_phi,
                       &psiw,
		       &psin,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  assert(ND(grad_phi) == 3);
  calculateGlobalExteriorNumericalFluxDarcySplitPressure(SHAPE(exteriorElementBoundaries)[0],
							 SHAPE(grad_phi)[1],
							 SHAPE(grad_phi)[2],
							 IDATA(exteriorElementBoundaries),
							 IDATA(elementBoundaryElements),
							 IDATA(elementBoundaryLocalElementBoundaries),
							 IDATA(isDOFBoundary),
							 DDATA(n),
							 DDATA(bc_a),
							 DDATA(bc_grad_phi),
							 DDATA(bc_psiw),
							 DDATA(bc_psin),
							 DDATA(a),
							 DDATA(grad_phi),
							 DDATA(psiw),
							 DDATA(psin),
							 DDATA(penalty),
							 DDATA(diffusiveFlux));
  
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcySplitPressure_sd(PyObject* self, 
									PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isDOFBoundary,*n,*bc_a,*bc_grad_phi,*bc_psiw,*bc_psin,
    *a,*grad_phi,*psiw,*psin,*penalty,*diffusiveFlux,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOO",
		       &rowptr,
		       &colind,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &isDOFBoundary,
                       &n,
                       &bc_a,
                       &bc_grad_phi,
                       &bc_psiw,
		       &bc_psin,
                       &a,
                       &grad_phi,
                       &psiw,
		       &psin,
                       &penalty,
                       &diffusiveFlux))
    return NULL;
  assert(ND(grad_phi) == 3);
  calculateGlobalExteriorNumericalFluxDarcySplitPressure_sd(SHAPE(exteriorElementBoundaries)[0],
							    SHAPE(grad_phi)[1],
							    SHAPE(grad_phi)[2],
							    IDATA(rowptr),
							    IDATA(colind),
							    IDATA(exteriorElementBoundaries),
							    IDATA(elementBoundaryElements),
							    IDATA(elementBoundaryLocalElementBoundaries),
							    IDATA(isDOFBoundary),
							    DDATA(n),
							    DDATA(bc_a),
							    DDATA(bc_grad_phi),
							    DDATA(bc_psiw),
							    DDATA(bc_psin),
							    DDATA(a),
							    DDATA(grad_phi),
							    DDATA(psiw),
							    DDATA(psin),
							    DDATA(penalty),
							    DDATA(diffusiveFlux));
 
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalFluxShallowWater_1D(PyObject* self,
                                                            PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*h,*hu,*advectiveFlux_h,*advectiveFlux_hu;
  double h_eps,tol_u,g;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOO",
                       &h_eps,
                       &tol_u,
		       &g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &h,
                       &hu,
		       &advectiveFlux_h,
		       &advectiveFlux_hu))
    return NULL;
  calculateInteriorNumericalFluxShallowWater_1D(SHAPE(interiorElementBoundaries)[0],
						SHAPE(hu)[1],
						SHAPE(hu)[2],
						h_eps,
                                                tol_u,
                                                g,
						IDATA(interiorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						DDATA(n),
						DDATA(h),
						DDATA(hu),
						DDATA(advectiveFlux_h),
						DDATA(advectiveFlux_hu));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateExteriorNumericalFluxShallowWater_1D(PyObject* self,
                                                            PyObject* args)
{
  PyObject *n,*h_l,*hu_l,*h_r,*hu_r,*advectiveFlux_h,*advectiveFlux_hu;
  double h_eps,tol_u,g;
  int nExteriorElementBoundaries;
  if(!PyArg_ParseTuple(args,"idddOOOOOOO",
		       &nExteriorElementBoundaries,
                       &h_eps,
                       &tol_u,
		       &g,
                       &n,
                       &h_l,
                       &hu_l,
                       &h_r,
                       &hu_r,
		       &advectiveFlux_h,
		       &advectiveFlux_hu))
    return NULL;
  calculateExteriorNumericalFluxShallowWater_1D(nExteriorElementBoundaries,
						SHAPE(hu_l)[1],
                                                h_eps,
                                                tol_u,
						g,
						DDATA(n),
						DDATA(h_l),
						DDATA(hu_l),
						DDATA(h_r),
						DDATA(hu_r),
						DDATA(advectiveFlux_h),
						DDATA(advectiveFlux_hu));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorNumericalFluxShallowWater_2D(PyObject* self,
								     PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*h,*hu,*hv,*advectiveFlux_h,*advectiveFlux_hu,*advectiveFlux_hv;
  double h_eps,tol_u,g;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOOOO",
                       &h_eps,
                       &tol_u,
		       &g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &h,
                       &hu,
                       &hv,
		       &advectiveFlux_h,
		       &advectiveFlux_hu,
		       &advectiveFlux_hv))
    return NULL;
  calculateInteriorNumericalFluxShallowWater_2D(SHAPE(interiorElementBoundaries)[0],
						SHAPE(hu)[1],
						SHAPE(hu)[2],
                                                h_eps,
                                                tol_u,
                                                g,
						IDATA(interiorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						DDATA(n),
						DDATA(h),
						DDATA(hu),
						DDATA(hv),
						DDATA(advectiveFlux_h),
						DDATA(advectiveFlux_hu),
						DDATA(advectiveFlux_hv));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalFluxShallowWater_2D(PyObject* self,
								     PyObject* args)
{
  PyObject *n,*h_l,*hu_l,*hv_l,*h_r,*hu_r,*hv_r,*advectiveFlux_h,*advectiveFlux_hu,*advectiveFlux_hv;
  double h_eps,tol_u,g;
  int nExteriorElementBoundaries;
  if(!PyArg_ParseTuple(args,"idddOOOOOOOOOO",
		       &nExteriorElementBoundaries,
                       &h_eps,
                       &tol_u,
		       &g,
                       &n,
                       &h_l,
                       &hu_l,
                       &hv_l,
                       &h_r,
                       &hu_r,
                       &hv_r,
		       &advectiveFlux_h,
		       &advectiveFlux_hu,
		       &advectiveFlux_hv))
    return NULL;
  calculateExteriorNumericalFluxShallowWater_2D(nExteriorElementBoundaries,
						SHAPE(hu_l)[1],
                                                h_eps,
                                                tol_u,
						g,
						DDATA(n),
						DDATA(h_l),
						DDATA(hu_l),
						DDATA(hv_l),
						DDATA(h_r),
						DDATA(hu_r),
						DDATA(hv_r),
						DDATA(advectiveFlux_h),
						DDATA(advectiveFlux_hu),
						DDATA(advectiveFlux_hv));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateInteriorNumericalFluxShallowWaterHLL_1D(PyObject* self,
							       PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*h,*hu,*advectiveFlux_h,*advectiveFlux_hu;
  double h_eps,tol_u,g;
  if(!PyArg_ParseTuple(args,"dddOOOOOOOO",
                       &h_eps,
                       &tol_u,
		       &g,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &h,
                       &hu,
		       &advectiveFlux_h,
		       &advectiveFlux_hu))
    return NULL;
  calculateInteriorNumericalFluxShallowWaterHLL_1D(SHAPE(interiorElementBoundaries)[0],
						   SHAPE(hu)[1],
						   SHAPE(hu)[2],
						   h_eps,
						   tol_u,
						   g,
						   IDATA(interiorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   DDATA(n),
						   DDATA(h),
						   DDATA(hu),
						   DDATA(advectiveFlux_h),
						   DDATA(advectiveFlux_hu));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cnumericalFluxCalculateExteriorNumericalFluxShallowWaterHLL_1D(PyObject* self,
							       PyObject* args)
{
  PyObject *n,*h_l,*hu_l,*h_r,*hu_r,*advectiveFlux_h,*advectiveFlux_hu;
  double h_eps,tol_u,g;
  int nExteriorElementBoundaries;
  if(!PyArg_ParseTuple(args,"idddOOOOOOO",
		       &nExteriorElementBoundaries,
                       &h_eps,
                       &tol_u,
		       &g,
                       &n,
                       &h_l,
                       &hu_l,
                       &h_r,
                       &hu_r,
		       &advectiveFlux_h,
		       &advectiveFlux_hu))
    return NULL;
  calculateExteriorNumericalFluxShallowWaterHLL_1D(nExteriorElementBoundaries,
						   SHAPE(hu_l)[1],
						   h_eps,
						   tol_u,
						   g,
						   DDATA(n),
						   DDATA(h_l),
						   DDATA(hu_l),
						   DDATA(h_r),
						   DDATA(hu_r),
						   DDATA(advectiveFlux_h),
						   DDATA(advectiveFlux_hu));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateInteriorChengShuNumericalFlux(PyObject* self,
						     PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *n,*u,*H,*dH,*H_element,*dH_element,
    *HamiltonJacobiFlux,*dHamiltonJacobiFlux_left,*dHamiltonJacobiFlux_right;
  int speedEvalFlag=0;
  if(!PyArg_ParseTuple(args,"iOOOOOOOOOOOO",
		       &speedEvalFlag,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &u,
                       &H,
                       &dH,
                       &H_element,
                       &dH_element,
                       &HamiltonJacobiFlux,
                       &dHamiltonJacobiFlux_left,
                       &dHamiltonJacobiFlux_right))
    return NULL;
  calculateInteriorChengShuNumericalFlux(SHAPE(interiorElementBoundaries)[0],
					 SHAPE(dH)[1],
					 SHAPE(dH)[2],
					 SHAPE(dH_element)[1],
					 SHAPE(dH)[3],
					 speedEvalFlag,
					 IDATA(interiorElementBoundaries),
					 IDATA(elementBoundaryElements),
					 IDATA(elementBoundaryLocalElementBoundaries),
					 DDATA(n),
					 DDATA(u),
					 DDATA(H),
					 DDATA(dH),
					 DDATA(H_element),
					 DDATA(dH_element),
					 DDATA(HamiltonJacobiFlux),
					 DDATA(dHamiltonJacobiFlux_left),
					 DDATA(dHamiltonJacobiFlux_right));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxApplySeepageFace(PyObject* self, 
			       PyObject* args)
{
  double eps;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isSeepageFace,*isDOFBoundary,*g,*n,*grad_u,*u,*advectiveFlux,*diffusiveFlux,*elementBoundaryDiameters;
  if(!PyArg_ParseTuple(args,"OOOOOdOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isSeepageFace,
		       &isDOFBoundary,
		       &eps,
		       &elementBoundaryDiameters,
		       &g,
		       &n,
		       &grad_u,
		       &u,
		       &advectiveFlux,
		       &diffusiveFlux))
    return NULL;
  applySeepageFace(SHAPE(exteriorElementBoundaries)[0],
		   SHAPE(grad_u)[1],
		   SHAPE(grad_u)[2],
		   IDATA(exteriorElementBoundaries),
		   IDATA(elementBoundaryElements),
		   IDATA(elementBoundaryLocalElementBoundaries),
		   IDATA(isSeepageFace),
		   IDATA(isDOFBoundary),
		   eps,
		   DDATA(elementBoundaryDiameters),
		   DDATA(g),
		   DDATA(n),
		   DDATA(grad_u),
		   DDATA(u),
		   DDATA(advectiveFlux),
		   DDATA(diffusiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxApplySeepageFaceJacobian(PyObject* self, 
				       PyObject* args)
{
  double eps;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*isSeepageFace,*g,*n,*grad_u,*u,*advectiveFlux,*diffusiveFlux,*v,*fluxJacobian,*elementBoundaryDiameters;
  if(!PyArg_ParseTuple(args,"OOOOdOOOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isSeepageFace,
		       &eps,
		       &elementBoundaryDiameters,
		       &g,
		       &n,
		       &grad_u,
		       &u,
		       &advectiveFlux,
		       &diffusiveFlux,
		       &v,
		       &fluxJacobian))
    return NULL;
  applySeepageFaceJacobian(SHAPE(exteriorElementBoundaries)[0],
			   SHAPE(grad_u)[1],
			   SHAPE(fluxJacobian)[2],
			   SHAPE(grad_u)[2],
			   IDATA(exteriorElementBoundaries),
			   IDATA(elementBoundaryElements),
			   IDATA(elementBoundaryLocalElementBoundaries),
			   IDATA(isSeepageFace),
			   eps,
			   DDATA(elementBoundaryDiameters),
			   DDATA(g),
			   DDATA(n),
			   DDATA(grad_u),
			   DDATA(u),
			   DDATA(advectiveFlux),
			   DDATA(diffusiveFlux),
			   DDATA(v),
			   DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateGlobalExteriorNumericalStressFlux(PyObject* self, 
							  PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *isDOFBoundary_w,
    *n,
    *bc_u,
    *bc_v,
    *bc_w,
    *stress,
    *u,
    *v,
    *w,
    *penalty,
    *stressFlux_u,
    *stressFlux_v,
    *stressFlux_w;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_u,
		       &isDOFBoundary_v,
		       &isDOFBoundary_w,
		       &n,
		       &bc_u,
		       &bc_v,
		       &bc_w,
		       &stress,
		       &u,
		       &v,
		       &w,
		       &penalty,
		       &stressFlux_u,
		       &stressFlux_v,
		       &stressFlux_w))
    return NULL;
  calculateGlobalExteriorNumericalStressFlux(SHAPE(n)[0],
					     SHAPE(n)[1],
					     SHAPE(n)[2],
					     IDATA(exteriorElementBoundaries),
					     IDATA(elementBoundaryElements),
					     IDATA(elementBoundaryLocalElementBoundaries),
					     IDATA(isDOFBoundary_u),
					     IDATA(isDOFBoundary_v),
					     IDATA(isDOFBoundary_w),
					     DDATA(n),
					     DDATA(bc_u),
					     DDATA(bc_v),
					     DDATA(bc_w),
					     DDATA(stress),
					     DDATA(u),
					     DDATA(v),
					     DDATA(w),
					     DDATA(penalty),
					     DDATA(stressFlux_u),
					     DDATA(stressFlux_v),
					     DDATA(stressFlux_w));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxUpdateExteriorNumericalStressFluxJacobian(PyObject* self, 
							PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *isDOFBoundary_u,
    *isDOFBoundary_v,
    *isDOFBoundary_w,
    *isStressBoundary_u,
    *isStressBoundary_v,
    *isStressBoundary_w,
    *n,
    *dstress_u_u,
    *dstress_u_v,
    *dstress_u_w,
    *dstress_v_u,
    *dstress_v_v,
    *dstress_v_w,
    *dstress_w_u,
    *dstress_w_v,
    *dstress_w_w,
    *v,
    *grad_v,
    *penalty,
    *fluxJacobian_u_u,
    *fluxJacobian_u_v,
    *fluxJacobian_u_w,
    *fluxJacobian_v_u,
    *fluxJacobian_v_v,
    *fluxJacobian_v_w,
    *fluxJacobian_w_u,
    *fluxJacobian_w_v,
    *fluxJacobian_w_w;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO",
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &isDOFBoundary_u,
		       &isDOFBoundary_v,
		       &isDOFBoundary_w,
		       &isStressBoundary_u,
		       &isStressBoundary_v,
		       &isStressBoundary_w,
		       &n,
		       &dstress_u_u,
		       &dstress_u_v,
		       &dstress_u_w,
		       &dstress_v_u,
		       &dstress_v_v,
		       &dstress_v_w,
		       &dstress_w_u,
		       &dstress_w_v,
		       &dstress_w_w,
		       &v,
		       &grad_v,
		       &penalty,
		       &fluxJacobian_u_u,
		       &fluxJacobian_u_v,
		       &fluxJacobian_u_w,
		       &fluxJacobian_v_u,
		       &fluxJacobian_v_v,
		       &fluxJacobian_v_w,
		       &fluxJacobian_w_u,
		       &fluxJacobian_w_v,
		       &fluxJacobian_w_w))
    return NULL;
  updateExteriorNumericalStressFluxJacobian(SHAPE(v)[0],
					    SHAPE(v)[1],
					    SHAPE(v)[2],
					    SHAPE(n)[2],
					    IDATA(exteriorElementBoundaries),
					    IDATA(elementBoundaryElements),
					    IDATA(elementBoundaryLocalElementBoundaries),
					    IDATA(isDOFBoundary_u),
					    IDATA(isDOFBoundary_v),
					    IDATA(isDOFBoundary_w),    
					    IDATA(isStressBoundary_u),
					    IDATA(isStressBoundary_v),
					    IDATA(isStressBoundary_w),
					    DDATA(n),
					    DDATA(dstress_u_u),
					    DDATA(dstress_u_v),
					    DDATA(dstress_u_w),
					    DDATA(dstress_v_u),
					    DDATA(dstress_v_v),
					    DDATA(dstress_v_w),
					    DDATA(dstress_w_u),
					    DDATA(dstress_w_v),
					    DDATA(dstress_w_w),
					    DDATA(v),
					    DDATA(grad_v),
					    DDATA(penalty),
					    DDATA(fluxJacobian_u_u),
					    DDATA(fluxJacobian_u_v),
					    DDATA(fluxJacobian_u_w),
					    DDATA(fluxJacobian_v_u),
					    DDATA(fluxJacobian_v_v),
					    DDATA(fluxJacobian_v_w),
					    DDATA(fluxJacobian_w_u),
					    DDATA(fluxJacobian_w_v),
					    DDATA(fluxJacobian_w_w));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalFluxRichards_sd(PyObject* self, 
							PyObject* args)
{
  PyObject *rowptr,
    *colind,
    *isSeepageFace,
    *isDOFBoundary,
    *n,
    *bc_u,
    *K,
    *grad_psi,
    *u,
    *K_rho_g,
    *penalty,
    *diffusiveFlux;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
                       &rowptr,
		       &colind,
		       &isSeepageFace,
		       &isDOFBoundary,
		       &n,
		       &bc_u,
		       &K,
		       &grad_psi,
		       &u,
		       &K_rho_g,
		       &penalty,
		       &diffusiveFlux))
    return NULL;
  calculateExteriorNumericalFluxRichards_sd(IDATA(rowptr),
					    IDATA(colind),
					    SHAPE(n)[0],
					    SHAPE(n)[1],
					    SHAPE(n)[2],
					    IDATA(isSeepageFace),
					    IDATA(isDOFBoundary),
					    DDATA(n),
					    DDATA(bc_u),
					    DDATA(K),
					    DDATA(grad_psi),
					    DDATA(u),
					    DDATA(K_rho_g),
					    DDATA(penalty),
					    DDATA(diffusiveFlux));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cnumericalFluxCalculateExteriorNumericalFluxJacobianRichards_sd(PyObject* self, 
								PyObject* args)
{
  PyObject *rowptr,
    *colind,
    *isDOFBoundary,
    *n,
    *bc_u,
    *K,
    *dK,
    *grad_psi,
    *grad_v,
    *u,
    *dK_rho_g,
    *v,
    *penalty,
    *fluxJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOO",
                       &rowptr,
		       &colind,
		       &isDOFBoundary,
		       &n,
		       &bc_u,
		       &K,
		       &dK,
		       &grad_psi,
		       &grad_v,
		       &u,
		       &dK_rho_g,
		       &v,
		       &penalty,
		       &fluxJacobian))
    return NULL;
  calculateExteriorNumericalFluxJacobianRichards_sd(IDATA(rowptr),
						    IDATA(colind),
						    SHAPE(grad_v)[0],
						    SHAPE(grad_v)[1],
						    SHAPE(grad_v)[2],
						    SHAPE(grad_v)[3],
						    IDATA(isDOFBoundary),
						    DDATA(n),
						    DDATA(bc_u),
						    DDATA(K),
						    DDATA(dK),
						    DDATA(grad_psi),
						    DDATA(grad_v),
						    DDATA(u),
						    DDATA(dK_rho_g),
						    DDATA(v),
						    DDATA(penalty),
						    DDATA(fluxJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef cnumericalFluxMethods[] = {
  {"calculateInteriorNumericalAdvectiveFluxConvexOneSonicPoint",
   cnumericalFluxCalculateInteriorNumericalAdvectiveFluxConvexOneSonicPoint,
   METH_VARARGS,
   "simple scalar riemann solve for convex fluxes and at most 1 sonic point"},
  {"calculateInteriorNumericalAdvectiveFluxRusanov",
   cnumericalFluxCalculateInteriorNumericalAdvectiveFluxRusanov,
   METH_VARARGS,
   "simple scalar riemann solve using Rusanov (local Lax-Friedrichs)"},
  {"calculateExteriorNumericalAdvectiveFluxRusanov",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxRusanov,
   METH_VARARGS,
   "simple scalar riemann solve using Rusanov (local Lax-Friedrichs) on boundaries"},
  {"calculateInteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound",
   cnumericalFluxCalculateInteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound,
   METH_VARARGS,
   "simple riemann solve using Rusanov (local Lax-Friedrichs) flux and eigenvalue bound that should be valid for a hyp. system"},
  {"calculateExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxRusanovWithEigenvalueBound,
   METH_VARARGS,
   "simple riemann solve using Rusanov (local Lax-Friedrichs) on boundaries and eigenvalue bound that should be valid for a hyp. system"},
  {"calculateExteriorNumericalAdvectiveFlux",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFlux,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries"},
  {"calculateExteriorNumericalAdvectiveFlux_free",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFlux_free,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries"},
  {"calculateExteriorNumericalAdvectiveFluxStokesP2D",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokesP2D,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries for Stokes equation"},
  {"calculateExteriorNumericalAdvectiveFluxStokesP3D",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokesP3D,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries for Stokes equation"},
  {"calculateExteriorNumericalAdvectiveFluxNavierStokes2D",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxNavierStokes2D,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries for NavierStokes equation"},
  {"calculateExteriorNumericalAdvectiveFluxNavierStokes3D",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxNavierStokes3D,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries for NavierStokes equation"},
  {"calculateExteriorNumericalAdvectiveFluxStokes2D",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokes2D,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries for Stokes equation"},
  {"calculateExteriorNumericalAdvectiveFluxStokes3D",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFluxStokes3D,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries for Stokes equation"},
  {"calculateInteriorNumericalAdvectiveFlux",
   cnumericalFluxCalculateInteriorNumericalAdvectiveFlux,
   METH_VARARGS,
   "calculate the numerical advective flux on interior element boundaries"},
  {"updateInteriorNumericalAdvectiveFluxJacobian",
   cnumericalFluxUpdateInteriorNumericalAdvectiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the advective flux on interior element boundaries"},
  {"updateInteriorTwoSidedNumericalFluxJacobian",
   cnumericalFluxUpdateInteriorTwoSidedNumericalFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of a two-sided numerical flux with the jacobian of the flux on interior element boundaries"},
  {"setInflowFlux",
   cnumericalFluxSetInflowFlux,
   METH_VARARGS,
   "set the inflow flux to the current flux"},
  {"calculateInteriorNumericalAdvectiveFlux_average",
   cnumericalFluxCalculateInteriorNumericalAdvectiveFlux_average,
   METH_VARARGS,
   "calculate the numerical advective flux on interior element boundaries"},
  {"calculateExteriorNumericalAdvectiveFlux_average",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFlux_average,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries"},
  {"calculateExteriorNumericalAdvectiveFlux_NoBC",
   cnumericalFluxCalculateExteriorNumericalAdvectiveFlux_NoBC,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries"},
  {"calculateGlobalExteriorInflowNumericalAdvectiveFlux",
   cnumericalFluxCalculateGlobalExteriorInflowNumericalAdvectiveFlux,
   METH_VARARGS,
   "calculate the numerical advective flux on global exterior element boundaries that are inflow"},
  {"updateExteriorNumericalAdvectiveFluxJacobian",
   cnumericalFluxUpdateExteriorNumericalAdvectiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the advective flux on exterior element boundaries"},
  {"updateExteriorNumericalAdvectiveFluxJacobian_free",
   cnumericalFluxUpdateExteriorNumericalAdvectiveFluxJacobian_free,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the advective flux on exterior element boundaries"},
  {"calculateInteriorNumericalDiffusiveFlux",
   cnumericalFluxCalculateInteriorNumericalDiffusiveFlux,
   METH_VARARGS,
   "calculate the numerical diffusive flux on interior element boundaries"},
  {"calculateInteriorNumericalDiffusiveFlux_sd",
   cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_sd,
   METH_VARARGS,
   "calculate the numerical diffusive flux on interior element boundaries"},
  {"updateInteriorNumericalDiffusiveFluxJacobian",
   cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the diffusive flux on interior element boundaries"},
  {"updateInteriorNumericalDiffusiveFluxJacobian_sd",
   cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_sd,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the jacobian of the diffusive flux on interior element boundaries"},
  {"calculateExteriorNumericalDiffusiveFlux",
   cnumericalFluxCalculateExteriorNumericalDiffusiveFlux,
   METH_VARARGS,
   "calculate the numerical diffusive flux on exterior element boundaries"},
  {"calculateExteriorNumericalDiffusiveFlux_sd",
   cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_sd,
   METH_VARARGS,
   "calculate the numerical diffusive flux on exterior element boundaries"},
  {"calculateExteriorNumericalDiffusiveFlux_free",
   cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_free,
   METH_VARARGS,
   "calculate the numerical diffusive flux on exterior element boundaries"},
  {"calculateExteriorNumericalDiffusiveFlux_free_sd",
   cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_free_sd,
   METH_VARARGS,
   "calculate the numerical diffusive flux on exterior element boundaries"},
  {"calculateExteriorNumericalDiffusiveFluxWithUpwinding_sd",
   cnumericalFluxCalculateExteriorNumericalDiffusiveFluxWithUpwinding_sd,
   METH_VARARGS,
   "calculate the numerical diffusive flux on exterior element boundaries but upwind based on potential"},
  {"updateExteriorNumericalDiffusiveFluxJacobian",
   cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the Jacobian of the diffusive flux on exterior element boundaries"},
  {"updateExteriorNumericalDiffusiveFluxJacobian_sd",
   cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_sd,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the Jacobian of the diffusive flux on exterior element boundaries"},
  {"updateExteriorNumericalDiffusiveFluxJacobian_free",
   cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_free,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the Jacobian of the diffusive flux on exterior element boundaries"},
  {"updateExteriorNumericalDiffusiveFluxJacobian_free_sd",
   cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_free_sd,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the Jacobian of the diffusive flux on exterior element boundaries"},
  {"updateExteriorNumericalDiffusiveFluxWithUpwindingJacobian_sd",
   cnumericalFluxUpdateExteriorNumericalDiffusiveFluxWithUpwindingJacobian_sd,
   METH_VARARGS,
   "update the Jacobian  of the numerical flux with the Jacobian of the diffusive flux on exterior element boundaries when upwinding diffusion"},
  { "calculateInteriorNumericalDiffusiveFlux_LDG_upwind",
    cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_LDG_upwind,
    METH_VARARGS, 
    "Calculate the diffusive flux on the interior element boundaries for LDG using and upwind approximation"},
  { "calculateInteriorNumericalDiffusiveFlux_LDG_upwind_sd",
    cnumericalFluxCalculateInteriorNumericalDiffusiveFlux_LDG_upwind_sd,
    METH_VARARGS, 
    "Calculate the diffusive flux on the interior element boundaries for LDG using and upwind approximation"},
  { "updateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind",
    cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind,
    METH_VARARGS, 
    "Update the diffusive flux Jacobian on the interior element boundaries for LDG using and upwind approximation"},
  { "updateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd",
    cnumericalFluxUpdateInteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd,
    METH_VARARGS, 
    "Update the diffusive flux Jacobian on the interior element boundaries for LDG using and upwind approximation"},
  { "calculateExteriorNumericalDiffusiveFlux_LDG_upwind",
    cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_LDG_upwind,
    METH_VARARGS, 
    "Calculate the diffusive flux on the exterior element boundaries for LDG using and upwind approximation"},
  { "calculateExteriorNumericalDiffusiveFlux_LDG_upwind_sd",
    cnumericalFluxCalculateExteriorNumericalDiffusiveFlux_LDG_upwind_sd,
    METH_VARARGS, 
    "Calculate the diffusive flux on the exterior element boundaries for LDG using and upwind approximation"},
  { "updateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind",
    cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind,
    METH_VARARGS, 
    "Update the diffusive flux Jacobian for LDG using and upwind approximation"},
  { "updateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd",
    cnumericalFluxUpdateExteriorNumericalDiffusiveFluxJacobian_LDG_upwind_sd,
    METH_VARARGS, 
    "Update the diffusive flux Jacobian for LDG using and upwind approximation"},
  { "calculateDiffusionMatrixSplittings_LDG_sd",
    cnumericalFluxCalculateDiffusionMatrixSplittings_LDG_sd,
    METH_VARARGS, 
    "Pick matrices in extended mixed formulation"},
  { "calculateInteriorLesaintRaviartNumericalFlux",
    cnumericalFluxCalculateInteriorLesaintRaviartNumericalFlux,
    METH_VARARGS, 
    "calculate  Lesaint Raviart boundary term as a nonconservative (2-sided) flux"},
  {"calculateExteriorLesaintRaviartNumericalFlux",
   cnumericalFluxCalculateExteriorLesaintRaviartNumericalFlux,
   METH_VARARGS,
   "calculate  Lesaint Raviart boundary term as a nonconservative (2-sided) flux on external boundaries"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCFF",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_fc_ff model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCFF_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_fc_ff model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for Twophase_fc_ff model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCFF_diffusiveFluxJacobian_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for Twophase_fc_ff model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFC",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_fc model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFC_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_fc model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for Twophase_fc model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFC_diffusiveFluxJacobian_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for Twophase_fc model"},
  {"calculateGlobalExteriorNumericalAdvectiveFlux_DarcyFC",
   cnumericalFluxCalculateGlobalExteriorNumericalAdvectiveFlux_DarcyFC,
   METH_VARARGS,
   "calculate the numerical advective flux on exterior element boundaries"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCPP",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_fc_pp model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCPP_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_fc_pp model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for Twophase_fc_pp model"},
  {"calculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcyFCPP_diffusiveFluxJacobian_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for Twophase_fc_pp model"},
  {"calculateGlobalExteriorNumericalFluxDarcySplitPressure",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcySplitPressure,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_split_pressure model"},
  {"calculateGlobalExteriorNumericalFluxDarcySplitPressure_sd",
   cnumericalFluxCalculateGlobalExteriorNumericalFluxDarcySplitPressure_sd,
   METH_VARARGS,
   "calculate weak dirichlet boundary terms for Twophase_split_pressure model"},
  {"calculateInteriorNumericalFluxShallowWater_1D",
   cnumericalFluxCalculateInteriorNumericalFluxShallowWater_1D,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for 1D Shallow Water Equations"},
  {"calculateExteriorNumericalFluxShallowWater_1D",
   cnumericalFluxCalculateExteriorNumericalFluxShallowWater_1D,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for 1D Shallow Water Equations"},
  {"calculateInteriorNumericalFluxShallowWater_2D",
   cnumericalFluxCalculateInteriorNumericalFluxShallowWater_2D,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for 2D Shallow Water Equations"},
  {"calculateExteriorNumericalFluxShallowWater_2D",
   cnumericalFluxCalculateExteriorNumericalFluxShallowWater_2D,
   METH_VARARGS,
   "calculate weak dirichlet boundary Jacobian terms for 2D Shallow Water Equations"},
  {"calculateInteriorNumericalFluxShallowWaterHLL_1D",
   cnumericalFluxCalculateInteriorNumericalFluxShallowWaterHLL_1D,
   METH_VARARGS,
   "HLL approximate Riemann solver for 1D Shallow Water Equations"},
  {"calculateExteriorNumericalFluxShallowWaterHLL_1D",
   cnumericalFluxCalculateExteriorNumericalFluxShallowWaterHLL_1D,
   METH_VARARGS,
   "HLL approximate Riemann solver for exterior boundaries and 1D Shallow Water Equations"},
  { "calculateInteriorChengShuNumericalFlux",
    cnumericalFluxCalculateInteriorChengShuNumericalFlux,
    METH_VARARGS, 
    "calculate Cheng-Shu flux term as a nonconservative (2-sided) flux"},
  { "applySeepageFace",
    cnumericalFluxApplySeepageFace,
    METH_VARARGS, 
    "apply a seepage face boundary condition"},
  { "applySeepageFaceJacobian",
    cnumericalFluxApplySeepageFaceJacobian,
    METH_VARARGS, 
    "apply a seepage face boundary condition (apply to jacobian)"},
  { "calculateGlobalExteriorNumericalStressFlux",
    cnumericalFluxCalculateGlobalExteriorNumericalStressFlux,
    METH_VARARGS, 
    "calculate the trace of the stress tensor for Dirichlet conditions"},
  { "updateExteriorNumericalStressFluxJacobian",
    cnumericalFluxUpdateExteriorNumericalStressFluxJacobian,
    METH_VARARGS, 
    "update the stress terms in the flux jacobians"},
  { "calculateExteriorNumericalFluxRichards_sd",
    cnumericalFluxCalculateExteriorNumericalFluxRichards_sd,
    METH_VARARGS, 
    "calculate the basic IIPG flux for Richards' equation"},
  { "calculateExteriorNumericalFluxJacobianRichards_sd",
    cnumericalFluxCalculateExteriorNumericalFluxJacobianRichards_sd,
    METH_VARARGS, 
    "calculate the basic IIPG flux Jacobian for Richards' equation"},
   { NULL,NULL,0,NULL}
};


PyMODINIT_FUNC initcnumericalFlux(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cnumericalFlux", cnumericalFluxMethods);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
