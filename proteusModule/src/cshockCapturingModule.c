#include "Python.h"
#include "numpy/arrayobject.h"
#include "shockCapturing.h"
/** \file cshockCapturingModule.c
    \defgroup cshockCapturing cshockCapturing
    \brief  Python interface to shockCapturing library
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))

static PyObject*
cshockCapturingCalculateNumericalDiffusionResGrad(PyObject* self,
                                             PyObject* args)
{
  double shockCapturingFactor;
  PyObject *elementDiameter,*strongResidual,*grad_u,*numDiff;
  if(!PyArg_ParseTuple(args,"dOOOO",
                       &shockCapturingFactor,
                       &elementDiameter,
                       &strongResidual,
                       &grad_u,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionResGrad(SHAPE(grad_u)[0],
                                     SHAPE(grad_u)[1],
                                     SHAPE(grad_u)[2],
                                     shockCapturingFactor,
                                     DDATA(elementDiameter),
                                     DDATA(strongResidual),
                                     DDATA(grad_u),
                                     DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cshockCapturingCalculateNumericalDiffusionResGradQuad(PyObject* self,
                                             PyObject* args)
{
  double shockCapturingFactor;
  PyObject *elementDiameter,*strongResidual,*grad_u,*numDiff;
  if(!PyArg_ParseTuple(args,"dOOOO",
                       &shockCapturingFactor,
                       &elementDiameter,
                       &strongResidual,
                       &grad_u,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionResGradQuad(SHAPE(grad_u)[0],
                                         SHAPE(grad_u)[1],
                                         SHAPE(grad_u)[2],
                                         shockCapturingFactor,
                                         DDATA(elementDiameter),
                                         DDATA(strongResidual),
                                         DDATA(grad_u),
                                         DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cshockCapturingCalculateNumericalDiffusionEikonal(PyObject* self,
						  PyObject* args)
{
  double shockCapturingFactor;
  PyObject *elementDiameter,*strongResidual,*numDiff;
  if(!PyArg_ParseTuple(args,"dOOO",
                       &shockCapturingFactor,
                       &elementDiameter,
                       &strongResidual,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionEikonal(SHAPE(numDiff)[0],
				     SHAPE(numDiff)[1],
				     shockCapturingFactor,
				     DDATA(elementDiameter),
				     DDATA(strongResidual),
				     DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cshockCapturingCalculateNumericalDiffusion_A_1(PyObject* self,
                                             PyObject* args)
{
  double shockCapturingFactor;
  PyObject *elementDiameter,*strongResidual,*mt,*df,*numDiff;
  if(!PyArg_ParseTuple(args,"dOOOOO",
                       &shockCapturingFactor,
                       &elementDiameter,
                       &strongResidual,
                       &mt,
                       &df,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusion_A_1(SHAPE(df)[0],
                                  SHAPE(df)[1],
                                  SHAPE(df)[2],
                                  shockCapturingFactor,
                                  DDATA(elementDiameter),
                                  DDATA(strongResidual),
                                  DDATA(mt),
                                  DDATA(df),
                                  DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cshockCapturingCalculateNumericalDiffusionHJ(PyObject* self,
					   PyObject* args)
{
  double shockCapturingFactor;
  char shockCapturingFlag;
  PyObject *elementDiameter,*strongResidual,*mt,*H,*numDiff;
  if(!PyArg_ParseTuple(args,"cdOOOOO",
		       &shockCapturingFlag,
                       &shockCapturingFactor,
                       &elementDiameter,
                       &strongResidual,
                       &mt,
                       &H,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionHJ(SHAPE(H)[0],
				SHAPE(H)[1],
				shockCapturingFlag,
				shockCapturingFactor,
				DDATA(elementDiameter),
				DDATA(strongResidual),
				DDATA(mt),
				DDATA(H),
				DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cshockCapturingCalculateNumericalDiffusionHJV2(PyObject* self,
					     PyObject* args)
{
  double shockCapturingFactor;
  char shockCapturingFlag;
  PyObject *elementDiameter,*strongResidual,*mt,*h,*numDiff;
  if(!PyArg_ParseTuple(args,"cdOOOOO",
		       &shockCapturingFlag,
                       &shockCapturingFactor,
                       &elementDiameter,
                       &strongResidual,
                       &mt,
                       &h,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionHJV2(SHAPE(h)[0],
				  SHAPE(h)[1],
				  shockCapturingFlag,
				  shockCapturingFactor,
				  DDATA(elementDiameter),
				  DDATA(strongResidual),
				  DDATA(mt),
				  DDATA(h),
				  DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cshockCapturingCalculateNumericalDiffusionResGradJuanes(PyObject* self,
							PyObject* args)
{
  double shockCapturingFactor,uSC;
  PyObject *elementDiameter,*strongResidual,*grad_u,*numDiff;
  if(!PyArg_ParseTuple(args,"ddOOOO",
                       &shockCapturingFactor,
		       &uSC,
                       &elementDiameter,
                       &strongResidual,
                       &grad_u,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionResGradJuanes(SHAPE(grad_u)[0],
					   SHAPE(grad_u)[1],
					   SHAPE(grad_u)[2],
					   shockCapturingFactor,
					   uSC,
					   DDATA(elementDiameter),
					   DDATA(strongResidual),
					   DDATA(grad_u),
					   DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cshockCapturingCalculateNumericalDiffusionJaffre(PyObject* self,
						 PyObject* args)
{
  double shockCapturingFactor,beta;
  PyObject *elementDiameter,*strongResidual,*grad_u,*numDiff;
  if(!PyArg_ParseTuple(args,"ddOOOO",
                       &shockCapturingFactor,
		       &beta,
                       &elementDiameter,
                       &strongResidual,
                       &grad_u,
                       &numDiff))
    return NULL;
  calculateNumericalDiffusionJaffre(SHAPE(grad_u)[0],
				    SHAPE(grad_u)[1],
				    SHAPE(grad_u)[2],
				    shockCapturingFactor,
				    beta,
				    DDATA(elementDiameter),
				    DDATA(strongResidual),
				    DDATA(grad_u),
				    DDATA(numDiff));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef cshockCapturingMethods[] = {
  { "calculateNumericalDiffusionResGrad", 
    cshockCapturingCalculateNumericalDiffusionResGrad,
    METH_VARARGS, 
    "Put in the standard shock capturing term for  linear  advection"},
  { "calculateNumericalDiffusionResGradQuad", 
    cshockCapturingCalculateNumericalDiffusionResGradQuad,
    METH_VARARGS, 
    "Put in the standard shock capturing term for  linear  advection"},
  { "calculateNumericalDiffusion_A_1", 
    cshockCapturingCalculateNumericalDiffusion_A_1,
    METH_VARARGS, 
    "Put in the standard shock capturing term for  linear  advection"},
  { "calculateNumericalDiffusionHJ", 
    cshockCapturingCalculateNumericalDiffusionHJ,
    METH_VARARGS, 
    "Put in shock capturing term for HJ from Barth and Sethian"},
  { "calculateNumericalDiffusionHJV2", 
    cshockCapturingCalculateNumericalDiffusionHJV2,
    METH_VARARGS, 
    "Put in shock capturing term for HJ play with formula some"},
  { "calculateNumericalDiffusionResGradJuanes", 
    cshockCapturingCalculateNumericalDiffusionResGradJuanes,
    METH_VARARGS, 
    "Put in Juanes version of shock capturing term"},
  { "calculateNumericalDiffusionJaffre", 
    cshockCapturingCalculateNumericalDiffusionJaffre,
    METH_VARARGS, 
    "Put in generic Jaffre formula for shock capturing term"},
  { "calculateNumericalDiffusionEikonal", 
    cshockCapturingCalculateNumericalDiffusionEikonal,
    METH_VARARGS, 
    "Put in a formula for shock capturing term for the Eikonal equation"},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcshockCapturing(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cshockCapturing", cshockCapturingMethods);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
