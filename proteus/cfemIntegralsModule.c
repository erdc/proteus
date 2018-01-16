#include "Python.h"
#include "numpy/arrayobject.h"
#include "femIntegrals.h"
#include "superluWrappersModule.h"
/** \defgroup cfemIntegralsModule cfemIntegralsModule
    \brief Python interface to femIntegrals library
    @{
*/
#define ND(p) ((PyArrayObject *)p)->nd
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define CSRVAL(p) ((double*)((SparseMatrix*)p)->A.nzval)

static PyObject* cfemIntegralsCopyLeftElementBoundaryInfo(PyObject* self,
                                                          PyObject* args)
{
  PyObject *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *exteriorElementBoundaries,
    *interiorElementBoundaries,
    *x,
    *n,
    *xg,
    *ng;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOO",
                        &elementBoundaryElements,
                        &elementBoundaryLocalElementBoundaries,
                        &exteriorElementBoundaries,
                        &interiorElementBoundaries,
                        &x,
                        &n,
                        &xg,
                        &ng))
    return NULL;
  copyLeftElementBoundaryInfo(SHAPE(n)[1],
                              SHAPE(n)[2],
                              SHAPE(n)[3],
                              SHAPE(exteriorElementBoundaries)[0],
                              SHAPE(interiorElementBoundaries)[0],
                              IDATA(elementBoundaryElements),
                              IDATA(elementBoundaryLocalElementBoundaries),
                              IDATA(exteriorElementBoundaries),
                              IDATA(interiorElementBoundaries),
                              DDATA(x),
                              DDATA(n),
                              DDATA(xg),
                              DDATA(ng));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsCopyGlobalElementBoundaryInfo(PyObject* self,
							    PyObject* args)
{
  PyObject *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *exteriorElementBoundaries,
    *interiorElementBoundaries,
    *x,
    *n,
    *ebqe_x,
    *ebqe_n,
    *xg,
    *ng;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOO",
                        &elementBoundaryElements,
                        &elementBoundaryLocalElementBoundaries,
                        &exteriorElementBoundaries,
                        &interiorElementBoundaries,
                        &x,
                        &n,
                        &ebqe_x,
                        &ebqe_n,
                        &xg,
                        &ng))
    return NULL;
  copyGlobalElementBoundaryInfo(SHAPE(n)[1],
				SHAPE(n)[2],
				SHAPE(n)[3],
				SHAPE(exteriorElementBoundaries)[0],
				SHAPE(interiorElementBoundaries)[0],
				IDATA(elementBoundaryElements),
				IDATA(elementBoundaryLocalElementBoundaries),
				IDATA(exteriorElementBoundaries),
				IDATA(interiorElementBoundaries),
				DDATA(x),
				DDATA(n),
				DDATA(ebqe_x),
				DDATA(ebqe_n),
				DDATA(xg),
				DDATA(ng));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsCopyLeftElementBoundaryInfo_movingDomain(PyObject* self,
								       PyObject* args)
{
  PyObject *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *exteriorElementBoundaries,
    *interiorElementBoundaries,
    *xt;
  if (!PyArg_ParseTuple(args,
                        "OOOOO",
                        &elementBoundaryElements,
                        &elementBoundaryLocalElementBoundaries,
                        &exteriorElementBoundaries,
                        &interiorElementBoundaries,
                        &xt))
    return NULL;
  copyLeftElementBoundaryInfo_movingDomain(SHAPE(xt)[1],
					   SHAPE(xt)[2],
					   SHAPE(exteriorElementBoundaries)[0],
					   SHAPE(interiorElementBoundaries)[0],
					   IDATA(elementBoundaryElements),
					   IDATA(elementBoundaryLocalElementBoundaries),
					   IDATA(exteriorElementBoundaries),
					   IDATA(interiorElementBoundaries),
					   DDATA(xt));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getValues(PyObject* self,
                                                                     PyObject* args)
{
  PyObject *psi,*v;
  if (!PyArg_ParseTuple(args,
                        "OO",
                        &psi,
                        &v))
    return NULL;
  parametricFiniteElementSpace_getValues(SHAPE(v)[0],
                                         SHAPE(v)[1],
                                         SHAPE(v)[2],
                                         DDATA(psi),
                                         DDATA(v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getValuesTrace(PyObject* self,
                                                                          PyObject* args)
{
  PyObject *psi,*permutations,*v;
  if (!PyArg_ParseTuple(args,
                        "OOO",
                        &psi,
                        &permutations,
                        &v))
    return NULL;
  parametricFiniteElementSpace_getValuesTrace(SHAPE(v)[0],
                                              SHAPE(v)[1],
                                              SHAPE(v)[2],
                                              SHAPE(v)[3],
                                              DDATA(psi),
                                              IDATA(permutations),
                                              DDATA(v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getValuesGlobalExteriorTrace(PyObject* self,
											PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *psi,*v;
  if (!PyArg_ParseTuple(args,
                        "OOOOO",
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
                        &psi,
                        &v))
    return NULL;
  parametricFiniteElementSpace_getValuesGlobalExteriorTrace(SHAPE(psi)[0],
							    SHAPE(psi)[1],
							    SHAPE(psi)[2],
							    SHAPE(exteriorElementBoundariesArray)[0],
							    IDATA(exteriorElementBoundariesArray),
							    IDATA(elementBoundaryElementsArray),
							    IDATA(elementBoundaryLocalElementBoundariesArray),
							    DDATA(psi),
							    DDATA(v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getGradientValues(PyObject* self,
                                                                             PyObject* args)
{
  PyObject *grad_psi,*inverseJacobian,*grad_v;
  if (!PyArg_ParseTuple(args,
                        "OOO",
                        &grad_psi,
                        &inverseJacobian,
                        &grad_v))
    return NULL;
  parametricFiniteElementSpace_getGradientValues(SHAPE(grad_v)[0],
                                                 SHAPE(grad_v)[1],
                                                 SHAPE(grad_v)[2],
                                                 SHAPE(grad_v)[3],
                                                 DDATA(grad_psi),
                                                 DDATA(inverseJacobian),
                                                 DDATA(grad_v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getHessianValues(PyObject* self,
									    PyObject* args)
{
  PyObject *Hessian_psi,*inverseJacobian,*Hessian_v;
  if (!PyArg_ParseTuple(args,
                        "OOO",
                        &Hessian_psi,
                        &inverseJacobian,
                        &Hessian_v))
    return NULL;
  parametricFiniteElementSpace_getHessianValues(SHAPE(Hessian_v)[0],
						SHAPE(Hessian_v)[1],
						SHAPE(Hessian_v)[2],
						SHAPE(Hessian_v)[3],
						DDATA(Hessian_psi),
						DDATA(inverseJacobian),
						DDATA(Hessian_v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getGradientValuesTrace(PyObject* self,
                                                                                  PyObject* args)
{
  PyObject *grad_psi,*permutations,*inverseJacobian,*grad_v;
  if (!PyArg_ParseTuple(args,
                        "OOOO",
                        &grad_psi,
                        &permutations,
                        &inverseJacobian,
                        &grad_v))
    return NULL;
  parametricFiniteElementSpace_getGradientValuesTrace(SHAPE(grad_v)[0],
                                                      SHAPE(grad_v)[1],
                                                      SHAPE(grad_v)[2],
                                                      SHAPE(grad_v)[3],
                                                      SHAPE(grad_v)[4],
                                                      DDATA(grad_psi),
                                                      IDATA(permutations),
                                                      DDATA(inverseJacobian),
                                                      DDATA(grad_v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricFiniteElementSpace_getGradientValuesGlobalExteriorTrace(PyObject* self,
												PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *grad_psi,*inverseJacobian,*grad_v;
  if (!PyArg_ParseTuple(args,
                        "OOOOOO",
                        &exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
			&grad_psi,
                        &inverseJacobian,
                        &grad_v))
    return NULL;
  parametricFiniteElementSpace_getGradientValuesGlobalExteriorTrace(SHAPE(grad_psi)[0],
								    SHAPE(grad_psi)[1],
								    SHAPE(grad_psi)[2],
								    SHAPE(grad_v)[3],
								    SHAPE(exteriorElementBoundariesArray)[0],
								    IDATA(exteriorElementBoundariesArray),
								    IDATA(elementBoundaryElementsArray),
								    IDATA(elementBoundaryLocalElementBoundariesArray),
								    DDATA(grad_psi),
								    DDATA(inverseJacobian),
								    DDATA(grad_v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getPermutations(PyObject* self,
                                                             PyObject* args)
{
  PyObject *xi,*permutations;
  if (!PyArg_ParseTuple(args,
                        "OO",
                        &xi,
                        &permutations))
    return NULL;
  parametricMaps_getPermutations(SHAPE(xi)[0],
                                 SHAPE(xi)[1],
                                 SHAPE(xi)[2],
                                 SHAPE(xi)[3],
                                 DDATA(xi),
                                 IDATA(permutations));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cfemIntegralsGetPermutationsGlobal(PyObject* self,
						    PyObject* args)
{
  PyObject *xArray,*xArrayNew,*permutations;
  if (!PyArg_ParseTuple(args,
                        "OOO",
                        &xArray,
			&xArrayNew,
                        &permutations))
    return NULL;
  getPermutationsGlobal(SHAPE(xArray)[0],
			SHAPE(xArray)[1],
			DDATA(xArray),
			DDATA(xArrayNew),
			IDATA(permutations));
  Py_INCREF(Py_None); 
  return Py_None;
}
/* static PyObject* cfemIntegralsParametricMaps_getPermutationsGlobalExterior(PyObject* self, */
/* 									   PyObject* args) */
/* { */
/*   PyObject *exteriorElementBoundariesArray, */
/*     *elementBoundaryElementsArray, */
/*     *elementBoundaryLocalElementBoundariesArray, */
/*     *xi,*permutations; */
/*   if (!PyArg_ParseTuple(args, */
/*                         "OOOOO", */
/* 			&exteriorElementBoundariesArray, */
/* 			&elementBoundaryElementsArray, */
/* 			&elementBoundaryLocalElementBoundariesArray, */
/*                         &xi, */
/*                         &permutations)) */
/*     return NULL; */
/*   parametricMaps_getPermutationsGlobalExterior(SHAPE(xi)[1], */
/* 					       SHAPE(xi)[2], */
/* 					       SHAPE(exteriorElementBoundariesArray)[0], */
/* 					       IDATA(exteriorElementBoundariesArray), */
/* 					       IDATA(elementBoundaryElementsArray), */
/* 					       IDATA(elementBoundaryLocalElementBoundariesArray), */
/* 					       DDATA(xi), */
/* 					       IDATA(permutations)); */
/*   Py_INCREF(Py_None);  */
/*   return Py_None; */
/* } */
static PyObject* cfemIntegralsParametricMaps_getValues(PyObject* self,
                                                       PyObject* args)
{
  PyObject *psi,*l2g,*nodes,*x;
  if (!PyArg_ParseTuple(args,
                        "OOOO",
                        &psi,
                        &l2g,
                        &nodes,
                        &x))
    return NULL;
  parametricMaps_getValues(SHAPE(x)[0],
                           SHAPE(x)[1],
                           SHAPE(l2g)[1],
                           SHAPE(x)[2],
                           DDATA(psi),
                           IDATA(l2g),
                           DDATA(nodes),
                           DDATA(x));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getValuesTrace(PyObject* self,
                                                            PyObject* args)
{
  PyObject *psi,*l2g,*nodes,*x;
  if (!PyArg_ParseTuple(args,
                        "OOOO",
                        &psi,
                        &l2g,
                        &nodes,
                        &x))
    return NULL;
  parametricMaps_getValuesTrace(SHAPE(x)[0],
                                SHAPE(x)[1],
                                SHAPE(x)[2],
                                SHAPE(l2g)[1],
                                SHAPE(x)[3],
                                DDATA(psi),
                                IDATA(l2g),
                                DDATA(nodes),
                                DDATA(x));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getValuesGlobalExteriorTrace(PyObject* self,
									  PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *psi,*l2g,
    *nodes,*x;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOO",
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
                        &psi,
                        &l2g,
                        &nodes,
                        &x))
    return NULL;
  parametricMaps_getValuesGlobalExteriorTrace(SHAPE(x)[1],
					      SHAPE(l2g)[1],
					      SHAPE(x)[2],
					      SHAPE(exteriorElementBoundariesArray)[0],
					      IDATA(exteriorElementBoundariesArray),
					      IDATA(elementBoundaryElementsArray),
					      IDATA(elementBoundaryLocalElementBoundariesArray),
					      DDATA(psi),
					      IDATA(l2g),
					      DDATA(nodes),
					      DDATA(x));
  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject* cfemIntegralsParametricMaps_getInverseValues(PyObject* self,
                                                              PyObject* args)
{
  PyObject *inverseJacobian,*l2g,*nodes,*x,*xi;
  if (!PyArg_ParseTuple(args,
                        "OOOOO",
                        &inverseJacobian,
                        &l2g,
                        &nodes,
                        &x,
                        &xi))
    return NULL;
  parametricMaps_getInverseValues(SHAPE(x)[0],
                                  SHAPE(x)[1],
                                  SHAPE(l2g)[1],
                                  SHAPE(inverseJacobian)[2],
                                  DDATA(inverseJacobian),
                                  IDATA(l2g),
                                  DDATA(nodes),
                                  DDATA(x),
                                  DDATA(xi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getInverseValuesTrace(PyObject* self,
                                                                   PyObject* args)
{
  PyObject *inverseJacobian,*l2g,*nodes,*x,*xi;
  if (!PyArg_ParseTuple(args,
                        "OOOOO",
                        &inverseJacobian,
                        &l2g,
                        &nodes,
                        &x,
                        &xi))
    return NULL;
  parametricMaps_getInverseValuesTrace(SHAPE(x)[0],
                                       SHAPE(x)[1],
                                       SHAPE(x)[2],
                                       SHAPE(l2g)[1],
                                       SHAPE(inverseJacobian)[3],
                                       DDATA(inverseJacobian),
                                       IDATA(l2g),
                                       DDATA(nodes),
                                       DDATA(x),
                                       DDATA(xi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getInverseValuesGlobalExteriorTrace(PyObject* self,
										 PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *inverseJacobian,*l2g,*nodes,*x,*xi;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOO",
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
                        &inverseJacobian,
                        &l2g,
                        &nodes,
                        &x,
                        &xi))
    return NULL;
  parametricMaps_getInverseValuesGlobalExteriorTrace(SHAPE(x)[1],
						     SHAPE(l2g)[1],
						     SHAPE(inverseJacobian)[3],
						     SHAPE(exteriorElementBoundariesArray)[0],
						     IDATA(exteriorElementBoundariesArray),
						     IDATA(elementBoundaryElementsArray),
						     IDATA(elementBoundaryLocalElementBoundariesArray),
						     DDATA(inverseJacobian),
						     IDATA(l2g),
						     DDATA(nodes),
						     DDATA(x),
						     DDATA(xi));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getJacobianValues(PyObject* self,
                                                               PyObject* args)
{
  PyObject *grad_psi,*l2g,*nodes,*jacobian,*jacobianDeterminant,*jacobianInverse;
  if (!PyArg_ParseTuple(args,
                        "OOOOOO",
                        &grad_psi,
                        &l2g,
                        &nodes,
                        &jacobian,
                        &jacobianDeterminant,
                        &jacobianInverse))
    return NULL;
  if (SHAPE(jacobian)[2] == 1)
    parametricMaps_getJacobianValues1D(SHAPE(jacobian)[0],
                                       SHAPE(jacobian)[1],
                                       SHAPE(l2g)[1],
                                       DDATA(grad_psi),
                                       IDATA(l2g),
                                       DDATA(nodes),
                                       DDATA(jacobian),
                                       DDATA(jacobianDeterminant),
                                       DDATA(jacobianInverse));
  else if (SHAPE(jacobian)[2] == 2)
    parametricMaps_getJacobianValues2D(SHAPE(jacobian)[0],
                                       SHAPE(jacobian)[1],
                                       SHAPE(l2g)[1],
                                       DDATA(grad_psi),
                                       IDATA(l2g),
                                       DDATA(nodes),
                                       DDATA(jacobian),
                                       DDATA(jacobianDeterminant),
                                       DDATA(jacobianInverse));
  else if (SHAPE(jacobian)[2] == 3)
    parametricMaps_getJacobianValues3D(SHAPE(jacobian)[0],
                                       SHAPE(jacobian)[1],
                                       SHAPE(l2g)[1],
                                       DDATA(grad_psi),
                                       IDATA(l2g),
                                       DDATA(nodes),
                                       DDATA(jacobian),
                                       DDATA(jacobianDeterminant),
                                       DDATA(jacobianInverse));
  else
    printf("error in getJacobianValues...jacobian not sized properly");
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getJacobianValuesTrace(PyObject* self,
                                                                    PyObject* args)
{
  PyObject *grad_psi,*boundaryNormals,*boundaryJacobians,*l2g,*nodes,*jacobianInverse,*metricTensor,*metricTensorDeterminantSqrt,*unitNormal;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOO",
                        &grad_psi,
                        &boundaryNormals,
                        &boundaryJacobians,
                        &l2g,
                        &nodes,
                        &jacobianInverse,
                        &metricTensor,
                        &metricTensorDeterminantSqrt,
                        &unitNormal))
    return NULL;
  if (SHAPE(jacobianInverse)[3] == 1)
    parametricMaps_getJacobianValuesTrace1D(SHAPE(jacobianInverse)[0],
                                            SHAPE(jacobianInverse)[1],
                                            SHAPE(jacobianInverse)[2],
                                            SHAPE(l2g)[1],
                                            DDATA(grad_psi),
                                            DDATA(boundaryNormals),
                                            DDATA(boundaryJacobians),
                                            IDATA(l2g),
                                            DDATA(nodes),
                                            DDATA(jacobianInverse),
                                            DDATA(metricTensor),
                                            DDATA(metricTensorDeterminantSqrt),
                                            DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[3] == 2)
    parametricMaps_getJacobianValuesTrace2D(SHAPE(jacobianInverse)[0],
                                            SHAPE(jacobianInverse)[1],
                                            SHAPE(jacobianInverse)[2],
                                            SHAPE(l2g)[1],
                                            DDATA(grad_psi),
                                            DDATA(boundaryNormals),
                                            DDATA(boundaryJacobians),
                                            IDATA(l2g),
                                            DDATA(nodes),
                                            DDATA(jacobianInverse),
                                            DDATA(metricTensor),
                                            DDATA(metricTensorDeterminantSqrt),
                                            DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[3] == 3)
    parametricMaps_getJacobianValuesTrace3D(SHAPE(jacobianInverse)[0],
                                            SHAPE(jacobianInverse)[1],
                                            SHAPE(jacobianInverse)[2],
                                            SHAPE(l2g)[1],
                                            DDATA(grad_psi),
                                            DDATA(boundaryNormals),
                                            DDATA(boundaryJacobians),
                                            IDATA(l2g),
                                            DDATA(nodes),
                                            DDATA(jacobianInverse),
                                            DDATA(metricTensor),
                                            DDATA(metricTensorDeterminantSqrt),
                                            DDATA(unitNormal));
  else
    printf("error in getJacobianValuesTrace...jacobianInverse not sized properly");
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getJacobianValuesTrace_movingDomain(PyObject* self,
										 PyObject* args)
{
  PyObject *xt,*grad_psi,*boundaryNormals,*boundaryJacobians,*l2g,*nodes,*jacobianInverse,*metricTensor,*metricTensorDeterminantSqrt,*unitNormal;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOO",
                        &xt,
			&grad_psi,
                        &boundaryNormals,
                        &boundaryJacobians,
                        &l2g,
                        &nodes,
                        &jacobianInverse,
                        &metricTensor,
                        &metricTensorDeterminantSqrt,
                        &unitNormal))
    return NULL;
  if (SHAPE(jacobianInverse)[3] == 1)
    parametricMaps_getJacobianValuesTrace1D(SHAPE(jacobianInverse)[0],
                                            SHAPE(jacobianInverse)[1],
                                            SHAPE(jacobianInverse)[2],
                                            SHAPE(l2g)[1],
                                            DDATA(grad_psi),
                                            DDATA(boundaryNormals),
                                            DDATA(boundaryJacobians),
                                            IDATA(l2g),
                                            DDATA(nodes),
                                            DDATA(jacobianInverse),
                                            DDATA(metricTensor),
                                            DDATA(metricTensorDeterminantSqrt),
                                            DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[3] == 2)
    parametricMaps_getJacobianValuesTrace2D_movingDomain(SHAPE(jacobianInverse)[0],
                                            SHAPE(jacobianInverse)[1],
                                            SHAPE(jacobianInverse)[2],
                                            SHAPE(l2g)[1],
							 DDATA(xt),
                                            DDATA(grad_psi),
                                            DDATA(boundaryNormals),
                                            DDATA(boundaryJacobians),
                                            IDATA(l2g),
                                            DDATA(nodes),
                                            DDATA(jacobianInverse),
                                            DDATA(metricTensor),
                                            DDATA(metricTensorDeterminantSqrt),
                                            DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[3] == 3)
    parametricMaps_getJacobianValuesTrace3D(SHAPE(jacobianInverse)[0],
                                            SHAPE(jacobianInverse)[1],
                                            SHAPE(jacobianInverse)[2],
                                            SHAPE(l2g)[1],
                                            DDATA(grad_psi),
                                            DDATA(boundaryNormals),
                                            DDATA(boundaryJacobians),
                                            IDATA(l2g),
                                            DDATA(nodes),
                                            DDATA(jacobianInverse),
                                            DDATA(metricTensor),
                                            DDATA(metricTensorDeterminantSqrt),
                                            DDATA(unitNormal));
  else
    printf("error in getJacobianValuesTrace...jacobianInverse not sized properly");
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getJacobianValuesGlobalExteriorTrace(PyObject* self,
										  PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *grad_psi,*boundaryNormals,*boundaryJacobians,*l2g,*nodes,*jacobianInverse,
    *metricTensor,*metricTensorDeterminantSqrt,*unitNormal;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOOO",
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
                        &grad_psi,
                        &boundaryNormals,
                        &boundaryJacobians,
                        &l2g,
                        &nodes,
                        &jacobianInverse,
                        &metricTensor,
                        &metricTensorDeterminantSqrt,
                        &unitNormal))
    return NULL;
  if (SHAPE(jacobianInverse)[2] == 1)
    parametricMaps_getJacobianValuesGlobalExteriorTrace1D(SHAPE(jacobianInverse)[1],
							  SHAPE(l2g)[1],
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
							  DDATA(grad_psi),
							  DDATA(boundaryNormals),
							  DDATA(boundaryJacobians),
							  IDATA(l2g),
							  DDATA(nodes),
							  DDATA(jacobianInverse),
							  DDATA(metricTensor),
							  DDATA(metricTensorDeterminantSqrt),
							  DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[2] == 2)
    parametricMaps_getJacobianValuesGlobalExteriorTrace2D(SHAPE(jacobianInverse)[1],
							  SHAPE(l2g)[1],
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
							  DDATA(grad_psi),
							  DDATA(boundaryNormals),
							  DDATA(boundaryJacobians),
							  IDATA(l2g),
							  DDATA(nodes),
							  DDATA(jacobianInverse),
							  DDATA(metricTensor),
							  DDATA(metricTensorDeterminantSqrt),
							  DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[2] == 3)
    parametricMaps_getJacobianValuesGlobalExteriorTrace3D(SHAPE(jacobianInverse)[1],
							  SHAPE(l2g)[1],
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
							  DDATA(grad_psi),
							  DDATA(boundaryNormals),
							  DDATA(boundaryJacobians),
							  IDATA(l2g),
							  DDATA(nodes),
							  DDATA(jacobianInverse),
							  DDATA(metricTensor),
							  DDATA(metricTensorDeterminantSqrt),
							  DDATA(unitNormal));
  else
    printf("error in getJacobianValuesTrace...jacobianInverse not sized properly");
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsParametricMaps_getJacobianValuesGlobalExteriorTrace_movingDomain(PyObject* self,
										  PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *xt,*grad_psi,*boundaryNormals,*boundaryJacobians,*l2g,*nodes,*jacobianInverse,
    *metricTensor,*metricTensorDeterminantSqrt,*unitNormal;
  if (!PyArg_ParseTuple(args,
                        "OOOOOOOOOOOOO",
			&exteriorElementBoundariesArray,
			&elementBoundaryElementsArray,
			&elementBoundaryLocalElementBoundariesArray,
                        &xt,
			&grad_psi,
                        &boundaryNormals,
                        &boundaryJacobians,
                        &l2g,
                        &nodes,
                        &jacobianInverse,
                        &metricTensor,
                        &metricTensorDeterminantSqrt,
                        &unitNormal))
    return NULL;
  if (SHAPE(jacobianInverse)[2] == 1)
    parametricMaps_getJacobianValuesGlobalExteriorTrace1D(SHAPE(jacobianInverse)[1],
							  SHAPE(l2g)[1],
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
							  DDATA(grad_psi),
							  DDATA(boundaryNormals),
							  DDATA(boundaryJacobians),
							  IDATA(l2g),
							  DDATA(nodes),
							  DDATA(jacobianInverse),
							  DDATA(metricTensor),
							  DDATA(metricTensorDeterminantSqrt),
							  DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[2] == 2)
    parametricMaps_getJacobianValuesGlobalExteriorTrace2D_movingDomain(SHAPE(jacobianInverse)[1],
							  SHAPE(l2g)[1],
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
								       DDATA(xt),
								       DDATA(grad_psi),
							  DDATA(boundaryNormals),
							  DDATA(boundaryJacobians),
							  IDATA(l2g),
							  DDATA(nodes),
							  DDATA(jacobianInverse),
							  DDATA(metricTensor),
							  DDATA(metricTensorDeterminantSqrt),
							  DDATA(unitNormal));
  else if (SHAPE(jacobianInverse)[2] == 3)
    parametricMaps_getJacobianValuesGlobalExteriorTrace3D(SHAPE(jacobianInverse)[1],
							  SHAPE(l2g)[1],
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
							  DDATA(grad_psi),
							  DDATA(boundaryNormals),
							  DDATA(boundaryJacobians),
							  IDATA(l2g),
							  DDATA(nodes),
							  DDATA(jacobianInverse),
							  DDATA(metricTensor),
							  DDATA(metricTensorDeterminantSqrt),
							  DDATA(unitNormal));
  else
    printf("error in getJacobianValuesTrace...jacobianInverse not sized properly");
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateMass_weak(PyObject* self, 
                                              PyObject* args)
{
  PyObject *mt,*w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &mt,
                       &w_dV,
                       &residual)) 
    return NULL;
  updateMass_weak(SHAPE(w_dV)[0],
                  SHAPE(w_dV)[1],
                  SHAPE(w_dV)[2],
                  DDATA(mt),
                  DDATA(w_dV),
                  DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateMassJacobian_weak(PyObject* self, 
                                                      PyObject* args)
{
  PyObject *dmt,*v_X_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dmt,
                       &v_X_w_dV,
                       &jacobian)) 
    return NULL;
  updateMassJacobian_weak(SHAPE(v_X_w_dV)[0],
                          SHAPE(v_X_w_dV)[1],
                          SHAPE(v_X_w_dV)[2],
                          SHAPE(v_X_w_dV)[3],
                          DDATA(dmt),
                          DDATA(v_X_w_dV),
                          DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateMassJacobian_weak_lowmem(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *dmt,*v,*w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dmt,
                       &v,
                       &w_dV,
                       &jacobian)) 
    return NULL;
  updateMassJacobian_weak_lowmem(SHAPE(v)[0],
                                 SHAPE(v)[1],
                                 SHAPE(v)[2],
                                 SHAPE(w_dV)[2],
                                 DDATA(dmt),
                                 DDATA(v),
                                 DDATA(w_dV),
                                 DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateMass_strong(PyObject* self, 
                                                PyObject* args)
{
  PyObject *mt,*strongResidual;
  if(!PyArg_ParseTuple(args,"OO",
                       &mt,
                       &strongResidual)) 
    return NULL;
  updateMass_strong(SHAPE(mt)[0],
                    SHAPE(mt)[1],
                    DDATA(mt),
                    DDATA(strongResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateMassJacobian_strong(PyObject* self, 
                                                        PyObject* args)
{
  PyObject *dmt,*v,*strongJacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dmt,
                       &v,
                       &strongJacobian)) 
    return NULL;
  updateMassJacobian_strong(SHAPE(v)[0],
                            SHAPE(v)[1],
                            SHAPE(v)[2],
                            DDATA(dmt),
                            DDATA(v),
                            DDATA(strongJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateMass_adjoint(PyObject* self, 
                                                 PyObject* args)
{
  PyObject *dmt,*w_dV,*Lstar_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dmt,
                       &w_dV,
                       &Lstar_w_dV)) 
    return NULL;
  updateMass_adjoint(SHAPE(w_dV)[0],
                     SHAPE(w_dV)[1],
                     SHAPE(w_dV)[2],
                     DDATA(dmt),
                     DDATA(w_dV),
                     DDATA(Lstar_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateAdvection_weak(PyObject* self, 
                                                   PyObject* args)
{
  PyObject *f,*w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &f,
                       &w_dV,
                       &residual))
    return NULL;
  updateAdvection_weak(SHAPE(w_dV)[0],
                       SHAPE(w_dV)[1],
                       SHAPE(w_dV)[2],
                       SHAPE(w_dV)[3],
                       DDATA(f),
                       DDATA(w_dV),
                       DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateAdvectionJacobian_weak(PyObject* self, 
                                                           PyObject* args)
{
  PyObject *df,*v_X_grad_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &df,
                       &v_X_grad_w_dV,
                       &jacobian))
    return NULL;
  updateAdvectionJacobian_weak(SHAPE(v_X_grad_w_dV)[0],
                               SHAPE(v_X_grad_w_dV)[1],
                               SHAPE(v_X_grad_w_dV)[2],
                               SHAPE(v_X_grad_w_dV)[3],
                               SHAPE(v_X_grad_w_dV)[4],
                               DDATA(df),
                               DDATA(v_X_grad_w_dV),
                               DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateAdvectionJacobian_weak_lowmem(PyObject* self, 
                                                                  PyObject* args)
{
  PyObject *df,*v,*grad_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &df,
                       &v,
                       &grad_w_dV,
                       &jacobian))
    return NULL;
  updateAdvectionJacobian_weak_lowmem(SHAPE(v)[0],
                                      SHAPE(v)[1],
                                      SHAPE(v)[2],
                                      SHAPE(grad_w_dV)[2],
                                      SHAPE(grad_w_dV)[3],
                                      DDATA(df),
                                      DDATA(v),
                                      DDATA(grad_w_dV),
                                      DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateAdvection_strong(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *df,*grad_u,*strongResidual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &df,
                       &grad_u,
                       &strongResidual))
    return NULL;
  updateAdvection_strong(SHAPE(grad_u)[0],
                         SHAPE(grad_u)[1],
                         SHAPE(grad_u)[2],
                         DDATA(df),
                         DDATA(grad_u),
                         DDATA(strongResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateAdvectionJacobian_strong(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *df,*grad_v,*strongJacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &df,
                       &grad_v,
                       &strongJacobian))
    return NULL;
  updateAdvectionJacobian_strong(SHAPE(grad_v)[0],
                                 SHAPE(grad_v)[1],
                                 SHAPE(grad_v)[2],
                                 SHAPE(grad_v)[3],
                                 DDATA(df),
                                 DDATA(grad_v),
                                 DDATA(strongJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateAdvection_adjoint(PyObject* self, 
                                                      PyObject* args)
{
  PyObject *df,*grad_w_dV,*Lstar_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &df,
                       &grad_w_dV,
                       &Lstar_w_dV))
    return NULL;
  updateAdvection_adjoint(SHAPE(grad_w_dV)[0],
                          SHAPE(grad_w_dV)[1],
                          SHAPE(grad_w_dV)[2],
                          SHAPE(grad_w_dV)[3],
                          DDATA(df),
                          DDATA(grad_w_dV),
                          DDATA(Lstar_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateHamiltonian_weak(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *H,*w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &H,
                       &w_dV,
                       &residual))
    return NULL;
  updateHamiltonian_weak(SHAPE(w_dV)[0],
                         SHAPE(w_dV)[1],
                         SHAPE(w_dV)[2],
                         DDATA(H),
                         DDATA(w_dV),
                         DDATA(residual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateHamiltonianJacobian_weak(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *dH,*grad_v_X_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dH,
                       &grad_v_X_w_dV,
                       &jacobian))
    return NULL;
  updateHamiltonianJacobian_weak(SHAPE(grad_v_X_w_dV)[0],
                                 SHAPE(grad_v_X_w_dV)[1],
                                 SHAPE(grad_v_X_w_dV)[2],
                                 SHAPE(grad_v_X_w_dV)[3],
                                 SHAPE(grad_v_X_w_dV)[4],
                                 DDATA(dH),
                                 DDATA(grad_v_X_w_dV),
                                 DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateHamiltonianJacobian_weak_lowmem(PyObject* self, 
                                                                    PyObject* args)
{
  PyObject *dH,*grad_v,*w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dH,
                       &grad_v,
                       &w_dV,
                       &jacobian))
    return NULL;
  updateHamiltonianJacobian_weak_lowmem(SHAPE(grad_v)[0],
                                        SHAPE(grad_v)[1],
                                        SHAPE(grad_v)[2],
                                        SHAPE(w_dV)[2],
                                        SHAPE(grad_v)[3],
                                        DDATA(dH),
                                        DDATA(grad_v),
                                        DDATA(w_dV),
                                        DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateHamiltonian_strong(PyObject* self, 
                                                       PyObject* args)
{
  PyObject *dH,*grad_u,*strongResidual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dH,
                       &grad_u,
                       &strongResidual))
    return NULL;
  updateHamiltonian_strong(SHAPE(grad_u)[0],
                           SHAPE(grad_u)[1],
                           SHAPE(grad_u)[2],
                           DDATA(dH),
                           DDATA(grad_u),
                           DDATA(strongResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateHamiltonianJacobian_strong(PyObject* self, 
                                                               PyObject* args)
{
  PyObject *dH,*grad_v,*strongJacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dH,
                       &grad_v,
                       &strongJacobian))
    return NULL;
  updateHamiltonianJacobian_strong(SHAPE(grad_v)[0],
                                   SHAPE(grad_v)[1],
                                   SHAPE(grad_v)[2],
                                   SHAPE(grad_v)[3],
                                   DDATA(dH),
                                   DDATA(grad_v),
                                   DDATA(strongJacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateHamiltonian_adjoint(PyObject* self, 
                                                        PyObject* args)
{
  PyObject *dH,*grad_w_dV,*Lstar_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dH,
                       &grad_w_dV,
                       &Lstar_w_dV))
    return NULL;
  updateHamiltonian_adjoint(SHAPE(grad_w_dV)[0],
                            SHAPE(grad_w_dV)[1],
                            SHAPE(grad_w_dV)[2],
                            SHAPE(grad_w_dV)[3],
                            DDATA(dH),
                            DDATA(grad_w_dV),
                            DDATA(Lstar_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_weak(PyObject* self, 
                                                   PyObject* args)
{
  PyObject *a,*grad_phi_X_grad_w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &a,
                       &grad_phi_X_grad_w_dV,
                       &residual))
    return NULL;
  updateDiffusion_weak(SHAPE(grad_phi_X_grad_w_dV)[0],
                       SHAPE(grad_phi_X_grad_w_dV)[1],
                       SHAPE(grad_phi_X_grad_w_dV)[2],
                       SHAPE(grad_phi_X_grad_w_dV)[3],
                       DDATA(a),
                       DDATA(grad_phi_X_grad_w_dV),
                       DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_weak_lowmem(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *a,*grad_phi,*grad_w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &a,
                       &grad_phi,
                       &grad_w_dV,
                       &residual))
    return NULL;
  updateDiffusion_weak_lowmem(SHAPE(grad_w_dV)[0],
                              SHAPE(grad_w_dV)[1],
                              SHAPE(grad_w_dV)[2],
                              SHAPE(grad_w_dV)[3],
                              DDATA(a),
                              DDATA(grad_phi),
                              DDATA(grad_w_dV),
                              DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_weak_sd(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *a,*grad_phi,*grad_w_dV,*residual,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &rowptr,
                       &colind,
                       &a,
                       &grad_phi,
                       &grad_w_dV,
                       &residual))
    return NULL;
  updateDiffusion_weak_sd(SHAPE(grad_w_dV)[0],
                          SHAPE(grad_w_dV)[1],
                          SHAPE(grad_w_dV)[2],
                          SHAPE(grad_w_dV)[3],
                          IDATA(rowptr),
                          IDATA(colind),
                          DDATA(a),
                          DDATA(grad_phi),
                          DDATA(grad_w_dV),
                          DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_weak(PyObject* self, 
                                                           PyObject* args)
{
  PyObject *l2g,*a,*da,*grad_phi_X_grad_w_dV,*dphi,*v,*grad_v_X_grad_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
                       &l2g,
                       &a,
                       &da,
                       &grad_phi_X_grad_w_dV,
                       &dphi,
                       &v,
                       &grad_v_X_grad_w_dV,
                       &jacobian))
    return NULL;
  updateDiffusionJacobian_weak(SHAPE(grad_v_X_grad_w_dV)[0],
                               SHAPE(grad_v_X_grad_w_dV)[1],
                               SHAPE(grad_v_X_grad_w_dV)[2],
                               SHAPE(grad_v_X_grad_w_dV)[3],
                               SHAPE(grad_v_X_grad_w_dV)[4],
                               IDATA(l2g),
                               DDATA(a),
                               DDATA(da),
                               DDATA(grad_phi_X_grad_w_dV),
                               DDATA(dphi),
                               DDATA(v),
                               DDATA(grad_v_X_grad_w_dV),
                               DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_weak_lowmem(PyObject* self, 
                                                                  PyObject* args)
{
  PyObject *l2g,*a,*da,*grad_phi,*grad_w_dV,*dphi,*v,*grad_v,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                       &l2g,
                       &a,
                       &da,
                       &grad_phi,
                       &grad_w_dV,
                       &dphi,
                       &v,
                       &grad_v,
                       &jacobian))
    return NULL;
  updateDiffusionJacobian_weak_lowmem(SHAPE(grad_w_dV)[0],
                                      SHAPE(grad_w_dV)[1],
                                      SHAPE(grad_v)[2],
                                      SHAPE(grad_w_dV)[2],
                                      SHAPE(grad_w_dV)[3],
                                      IDATA(l2g),
                                      DDATA(a),
                                      DDATA(da),
                                      DDATA(grad_phi),
                                      DDATA(grad_w_dV),
                                      DDATA(dphi),
                                      DDATA(v),
                                      DDATA(grad_v),
                                      DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_weak_sd(PyObject* self, 
                                                                  PyObject* args)
{
  PyObject *l2g,*a,*da,*grad_phi,*grad_w_dV,*dphi,*v,*grad_v,*jacobian,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
                       &rowptr,
                       &colind,
                       &l2g,
                       &a,
                       &da,
                       &grad_phi,
                       &grad_w_dV,
                       &dphi,
                       &v,
                       &grad_v,
                       &jacobian))
    return NULL;
  updateDiffusionJacobian_weak_sd(SHAPE(grad_w_dV)[0],
                                  SHAPE(grad_w_dV)[1],
                                  SHAPE(grad_v)[2],
                                  SHAPE(grad_w_dV)[2],
                                  SHAPE(grad_w_dV)[3],
                                  IDATA(rowptr),
                                  IDATA(colind),
                                  IDATA(l2g),
                                  DDATA(a),
                                  DDATA(da),
                                  DDATA(grad_phi),
                                  DDATA(grad_w_dV),
                                  DDATA(dphi),
                                  DDATA(v),
                                  DDATA(grad_v),
                                  DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_strong(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *da,*grad_phi,*grad_u,*strongResidual;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &da,
                       &grad_phi,
                       &grad_u,
                       &strongResidual))
    return NULL;
  updateDiffusion_strong(SHAPE(grad_u)[0],
                         SHAPE(grad_u)[1],
                         SHAPE(grad_u)[2],
                         DDATA(da),
                         DDATA(grad_phi),
                         DDATA(grad_u),
                         DDATA(strongResidual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_strong_sd(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *da,*grad_phi,*grad_u,*strongResidual,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &rowptr,
                       &colind,
                       &da,
                       &grad_phi,
                       &grad_u,
                       &strongResidual))
    return NULL;
  updateDiffusion_strong_sd(SHAPE(grad_u)[0],
                            SHAPE(grad_u)[1],
                            SHAPE(grad_u)[2],
                            IDATA(rowptr),
                            IDATA(colind),
                            DDATA(da),
                            DDATA(grad_phi),
                            DDATA(grad_u),
                            DDATA(strongResidual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_strong(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *l2g,*da,*dphi,*grad_phi,*grad_u,*grad_v,*strongJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
                       &l2g,
                       &da,
                       &dphi,
                       &grad_phi,
                       &grad_u,
                       &grad_v,
                       &strongJacobian))
    return NULL;
  updateDiffusionJacobian_strong(SHAPE(grad_v)[0],
                                 SHAPE(grad_v)[1],
                                 SHAPE(grad_v)[2],
                                 SHAPE(grad_v)[3],
                                 IDATA(l2g),
                                 DDATA(da),
                                 DDATA(dphi),
                                 DDATA(grad_phi),
                                 DDATA(grad_u),
                                 DDATA(grad_v),
                                 DDATA(strongJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_strong_sd(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *l2g,*da,*dphi,*grad_phi,*grad_u,*grad_v,*strongJacobian,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                       &rowptr,
                       &colind,
                       &l2g,
                       &da,
                       &dphi,
                       &grad_phi,
                       &grad_u,
                       &grad_v,
                       &strongJacobian))
    return NULL;
  updateDiffusionJacobian_strong_sd(SHAPE(grad_v)[0],
                                    SHAPE(grad_v)[1],
                                    SHAPE(grad_v)[2],
                                    SHAPE(grad_v)[3],
                                    IDATA(rowptr),
                                    IDATA(colind),
                                    IDATA(l2g),
                                    DDATA(da),
                                    DDATA(dphi),
                                    DDATA(grad_phi),
                                    DDATA(grad_u),
                                    DDATA(grad_v),
                                    DDATA(strongJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion2_strong(PyObject* self, 
						      PyObject* args)
{
  PyObject *a,*Hess_phi,*strongResidual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &a,
                       &Hess_phi,
                       &strongResidual))
    return NULL;
  updateDiffusion2_strong(SHAPE(Hess_phi)[0],
			  SHAPE(Hess_phi)[1],
			  SHAPE(Hess_phi)[2],
			  DDATA(a),
			  DDATA(Hess_phi),
			  DDATA(strongResidual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion2_strong_sd(PyObject* self, 
                                                         PyObject* args)
{
  PyObject *a,*Hess_phi,*strongResidual,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &rowptr,
                       &colind,
                       &a,
                       &Hess_phi,
                       &strongResidual))
    return NULL;
  updateDiffusion2_strong_sd(SHAPE(Hess_phi)[0],
                             SHAPE(Hess_phi)[1],
                             SHAPE(Hess_phi)[2],
                             IDATA(rowptr),
                             IDATA(colind),
                             DDATA(a),
                             DDATA(Hess_phi),
                             DDATA(strongResidual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian2_strong(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *l2g,*a,*da,*dphi,*Hess_phi,*v,*Hess_v,*strongJacobian;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
                       &l2g,
		       &a,
                       &da,
		       &v,
                       &Hess_phi,
                       &dphi,
                       &Hess_v,
                       &strongJacobian))
    return NULL;
  updateDiffusionJacobian2_strong(SHAPE(Hess_v)[0],
				  SHAPE(Hess_v)[1],
				  SHAPE(Hess_v)[2],
				  SHAPE(Hess_v)[3],
				  IDATA(l2g),
				  DDATA(a),
				  DDATA(da),
				  DDATA(v),
				  DDATA(Hess_phi),
				  DDATA(dphi),
				  DDATA(Hess_v),
				  DDATA(strongJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian2_strong_sd(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *l2g,*a,*da,*dphi,*Hess_phi,*v,*Hess_v,*strongJacobian,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
                       &rowptr,
                       &colind,
                       &l2g,
		       &a,
                       &da,
		       &v,
                       &Hess_phi,
                       &dphi,
                       &Hess_v,
                       &strongJacobian))
    return NULL;
  updateDiffusionJacobian2_strong_sd(SHAPE(Hess_v)[0],
                                     SHAPE(Hess_v)[1],
                                     SHAPE(Hess_v)[2],
                                     SHAPE(Hess_v)[3],
                                     IDATA(rowptr),
                                     IDATA(colind),
                                     IDATA(l2g),
                                     DDATA(a),
                                     DDATA(da),
                                     DDATA(v),
                                     DDATA(Hess_phi),
                                     DDATA(dphi),
                                     DDATA(Hess_v),
                                     DDATA(strongJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_adjoint(PyObject* self, 
                                                      PyObject* args)
{
  PyObject *da,*grad_phi,*grad_w_dV,*Lstar_w_dV;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &da,
                       &grad_phi,
                       &grad_w_dV,
                       &Lstar_w_dV))
    return NULL;
  updateDiffusion_adjoint(SHAPE(grad_w_dV)[0],
                          SHAPE(grad_w_dV)[1],
                          SHAPE(grad_w_dV)[2],
                          SHAPE(grad_w_dV)[3],
                          DDATA(da),
                          DDATA(grad_phi),
                          DDATA(grad_w_dV),
                          DDATA(Lstar_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_adjoint_sd(PyObject* self, 
                                                         PyObject* args)
{
  PyObject *da,*grad_phi,*grad_w_dV,*Lstar_w_dV,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &rowptr,
                       &colind,
                       &da,
                       &grad_phi,
                       &grad_w_dV,
                       &Lstar_w_dV))
    return NULL;
  updateDiffusion_adjoint_sd(SHAPE(grad_w_dV)[0],
                             SHAPE(grad_w_dV)[1],
                             SHAPE(grad_w_dV)[2],
                             SHAPE(grad_w_dV)[3],
                             IDATA(rowptr),
                             IDATA(colind),
                             DDATA(da),
                             DDATA(grad_phi),
                             DDATA(grad_w_dV),
                             DDATA(Lstar_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion2_adjoint(PyObject* self, 
						       PyObject* args)
{
  PyObject *a,*Hess_w_dV,*Lstar_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &a,
                       &Hess_w_dV,
                       &Lstar_w_dV))
    return NULL;
  updateDiffusion2_adjoint(SHAPE(Hess_w_dV)[0],
			   SHAPE(Hess_w_dV)[1],
			   SHAPE(Hess_w_dV)[2],
			   SHAPE(Hess_w_dV)[3],
			   DDATA(a),
			   DDATA(Hess_w_dV),
			   DDATA(Lstar_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion2_adjoint_sd(PyObject* self, 
                                                          PyObject* args)
{
  PyObject *a,*Hess_w_dV,*Lstar_w_dV,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &rowptr,
                       &colind,
                       &a,
                       &Hess_w_dV,
                       &Lstar_w_dV))
    return NULL;
  updateDiffusion2_adjoint_sd(SHAPE(Hess_w_dV)[0],
                              SHAPE(Hess_w_dV)[1],
                              SHAPE(Hess_w_dV)[2],
                              SHAPE(Hess_w_dV)[3],
                              IDATA(rowptr),
                              IDATA(colind),
                              DDATA(a),
                              DDATA(Hess_w_dV),
                              DDATA(Lstar_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateReaction_weak(PyObject* self, 
                                 PyObject* args)
{
  PyObject *r,*w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &r,
                       &w_dV,
                       &residual)) 
    return NULL;
  updateReaction_weak(SHAPE(w_dV)[0],
                      SHAPE(w_dV)[1],
                      SHAPE(w_dV)[2],
                      DDATA(r),
                      DDATA(w_dV),
                      DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
cfemIntegralsUpdateReactionJacobian_weak(PyObject* self, 
                                         PyObject* args)
{
  PyObject *dr,*v_X_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dr,
                       &v_X_w_dV,
                       &jacobian)) 
    return NULL;
  updateReactionJacobian_weak(SHAPE(v_X_w_dV)[0],
                              SHAPE(v_X_w_dV)[1],
                              SHAPE(v_X_w_dV)[2],
                              SHAPE(v_X_w_dV)[3],
                              DDATA(dr),
                              DDATA(v_X_w_dV),
                              DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
cfemIntegralsUpdateReactionJacobian_weak_lowmem(PyObject* self, 
                                                PyObject* args)
{
  PyObject *dr,*v,*w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dr,
                       &v,
                       &w_dV,
                       &jacobian)) 
    return NULL;
  updateReactionJacobian_weak_lowmem(SHAPE(v)[0],
                                     SHAPE(v)[1],
                                     SHAPE(v)[2],
                                     SHAPE(w_dV)[2],
                                     DDATA(dr),
                                     DDATA(v),
                                     DDATA(w_dV),
                                     DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateReaction_strong(PyObject* self, 
                                   PyObject* args)
{
  PyObject *r,*strongResidual;
  if(!PyArg_ParseTuple(args,"OO",
                       &r,
                       &strongResidual)) 
    return NULL;
  updateReaction_strong(SHAPE(r)[0],
                        SHAPE(r)[1],
                        DDATA(r),
                        DDATA(strongResidual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* 
cfemIntegralsUpdateReactionJacobian_strong(PyObject* self, 
                                           PyObject* args)
{
  PyObject *dr,*v,*strongJacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dr,
                       &v,
                       &strongJacobian)) 
    return NULL;
  updateReactionJacobian_strong(SHAPE(v)[0],
                                SHAPE(v)[1],
                                SHAPE(v)[2],
                                DDATA(dr),
                                DDATA(v),
                                DDATA(strongJacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateReaction_adjoint(PyObject* self, 
                                    PyObject* args)
{
  PyObject *dr,*w_dV,*Lstar_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &dr,
                       &w_dV,
                       &Lstar_w_dV)) 
    return NULL;
  updateReaction_adjoint(SHAPE(w_dV)[0],
                         SHAPE(w_dV)[1],
                         SHAPE(w_dV)[2],
                         DDATA(dr),
                         DDATA(w_dV),
                         DDATA(Lstar_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateSubgridError(PyObject* self, 
                                PyObject* args)
{
  PyObject *error,*Lstar_w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &error,
                       &Lstar_w_dV,
                       &residual)) 
    return NULL;
  updateSubgridError(SHAPE(Lstar_w_dV)[0],
                     SHAPE(Lstar_w_dV)[1],
                     SHAPE(Lstar_w_dV)[2],
                     DDATA(error),
                     DDATA(Lstar_w_dV),
                     DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateSubgridErrorJacobian(PyObject* self, 
                                        PyObject* args)
{
  PyObject *derror,*Lstar_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &derror,
                       &Lstar_w_dV,
                       &jacobian)) 
    return NULL;
  updateSubgridErrorJacobian(SHAPE(Lstar_w_dV)[0],
                             SHAPE(Lstar_w_dV)[1],
                             SHAPE(derror)[2],
                             SHAPE(Lstar_w_dV)[2],
                             DDATA(derror),
                             DDATA(Lstar_w_dV),
                             DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateNumericalDiffusion(PyObject* self, 
                                      PyObject* args)
{
  PyObject *numDiff,*grad_u_X_grad_w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOO",
                       &numDiff,
                       &grad_u_X_grad_w_dV,
                       &residual)) 
    return NULL;
  updateNumericalDiffusion(SHAPE(grad_u_X_grad_w_dV)[0],
                           SHAPE(grad_u_X_grad_w_dV)[1],
                           SHAPE(grad_u_X_grad_w_dV)[2],
                           SHAPE(grad_u_X_grad_w_dV)[3],
                           DDATA(numDiff),
                           DDATA(grad_u_X_grad_w_dV),
                           DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateNumericalDiffusion_lowmem(PyObject* self, 
                                      PyObject* args)
{
  PyObject *numDiff,*grad_u,*grad_w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &numDiff,
                       &grad_u,
                       &grad_w_dV,
                       &residual)) 
    return NULL;
  updateNumericalDiffusion_lowmem(SHAPE(grad_w_dV)[0],
                                  SHAPE(grad_w_dV)[1],
                                  SHAPE(grad_w_dV)[2],
                                  SHAPE(grad_w_dV)[3],
                                  DDATA(numDiff),
                                  DDATA(grad_u),
                                  DDATA(grad_w_dV),
                                  DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateNumericalDiffusionJacobian(PyObject* self, 
                                              PyObject* args)
{
  PyObject *numDiff,*grad_v_X_grad_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOO",
                       &numDiff,
                       &grad_v_X_grad_w_dV,
                       &jacobian)) 
    return NULL;
  updateNumericalDiffusionJacobian(SHAPE(grad_v_X_grad_w_dV)[0],
                                   SHAPE(grad_v_X_grad_w_dV)[1],
                                   SHAPE(grad_v_X_grad_w_dV)[2],
                                   SHAPE(grad_v_X_grad_w_dV)[3],
                                   SHAPE(grad_v_X_grad_w_dV)[4],
                                   DDATA(numDiff),
                                   DDATA(grad_v_X_grad_w_dV),
                                   DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateNumericalDiffusionJacobian_lowmem(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *numDiff,*grad_v,*grad_w_dV,*jacobian;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &numDiff,
                       &grad_v,
                       &grad_w_dV,
                       &jacobian)) 
    return NULL;
  updateNumericalDiffusionJacobian_lowmem(SHAPE(grad_v)[0],
                                          SHAPE(grad_v)[1],
                                          SHAPE(grad_v)[2],
                                          SHAPE(grad_w_dV)[2],
                                          SHAPE(grad_w_dV)[3],
                                          DDATA(numDiff),
                                          DDATA(grad_v),
                                          DDATA(grad_w_dV),
                                          DDATA(jacobian));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateScalarScalarProduct(PyObject* self, 
                                          PyObject* args)
{
  PyObject *sL,*sR,*sResult;
  if(!PyArg_ParseTuple(args,"OOO",
                       &sL,
                       &sR,
                       &sResult))
    return NULL;
  calculateScalarScalarProduct(SHAPE(sL)[0],
                               SHAPE(sL)[1],
                               DDATA(sL),
                               DDATA(sR),
                               DDATA(sResult));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateVectorScalarProduct(PyObject* self, 
                                          PyObject* args)
{
  PyObject *v,*s,*vResult;
  if(!PyArg_ParseTuple(args,"OOO",
                       &v,
                       &s,
                       &vResult)) 
    return NULL;
  calculateVectorScalarProduct(SHAPE(v)[0],
                               SHAPE(v)[1],
                               SHAPE(v)[2],
                               DDATA(v),
                               DDATA(s),
                               DDATA(vResult));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateTensorScalarProduct(PyObject* self, 
                                          PyObject* args)
{
  PyObject *t,*s,*tResult;
  if(!PyArg_ParseTuple(args,"OOO",
                       &t,
                       &s,
                       &tResult)) 
    return NULL;
  calculateTensorScalarProduct(SHAPE(t)[0],
                               SHAPE(t)[1],
                               SHAPE(t)[2],
                               DDATA(t),
                               DDATA(s),
                               DDATA(tResult));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedShape(PyObject* self,
                                    PyObject* args)
{
  PyObject *dVR,*abs_det_jac,*w,*w_dV;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dVR,
                       &abs_det_jac,
                       &w,
                       &w_dV))
    return NULL;
  calculateWeightedShape(SHAPE(w_dV)[0], /* nElements_global */
                         SHAPE(w_dV)[1], /* nQuadraturePoints_element*/
                         SHAPE(w_dV)[2], /* nDOF_test_element */
                         DDATA(dVR),
                         DDATA(abs_det_jac),
                         DDATA(w),
                         DDATA(w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedShapeGradients(PyObject* self,
                                             PyObject* args)
{
  PyObject *dVR,*abs_det_jac,*grad_w,*grad_w_dV;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dVR,
                       &abs_det_jac,
                       &grad_w,
                       &grad_w_dV))
    return NULL;
  calculateWeightedShapeGradients(SHAPE(grad_w_dV)[0],
                                  SHAPE(grad_w_dV)[1],
                                  SHAPE(grad_w_dV)[2],
                                  SHAPE(grad_w_dV)[3],
                                  DDATA(dVR),
                                  DDATA(abs_det_jac),
                                  DDATA(grad_w),
                                  DDATA(grad_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedPiolaShapeGradients(PyObject* self,
						  PyObject* args)
{
  PyObject *dVR,*abs_det_jac,*grad_w,*grad_w_dV;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dVR,
                       &abs_det_jac,
                       &grad_w,
                       &grad_w_dV))
    return NULL;
  calculateWeightedPiolaShapeGradients(SHAPE(grad_w_dV)[0],
				       SHAPE(grad_w_dV)[1],
				       SHAPE(grad_w_dV)[2],
				       SHAPE(grad_w_dV)[3],
				       DDATA(dVR),
				       DDATA(abs_det_jac),
				       DDATA(grad_w),
				       DDATA(grad_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedShapeHessians(PyObject* self,
                                             PyObject* args)
{
  PyObject *dVR,*abs_det_jac,*Hess_w,*Hess_w_dV;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dVR,
                       &abs_det_jac,
                       &Hess_w,
                       &Hess_w_dV))
    return NULL;
  calculateWeightedShapeHessians(SHAPE(Hess_w_dV)[0],
				 SHAPE(Hess_w_dV)[1],
				 SHAPE(Hess_w_dV)[2],
				 SHAPE(Hess_w_dV)[3],
				 DDATA(dVR),
				 DDATA(abs_det_jac),
				 DDATA(Hess_w),
				 DDATA(Hess_w_dV));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateShape_X_weightedShape(PyObject* self, 
                                            PyObject* args)
{
  PyObject *v,*w_dV,*v_X_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &v,
                       &w_dV,
                       &v_X_w_dV))
    return NULL;
  calculateShape_X_weightedShape(SHAPE(v_X_w_dV)[0],
                                 SHAPE(v_X_w_dV)[1],
                                 SHAPE(v_X_w_dV)[2],
                                 SHAPE(v_X_w_dV)[3],
                                 DDATA(v),
                                 DDATA(w_dV),
                                 DDATA(v_X_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateShape_X_weightedGradShape(PyObject* self, 
                                                PyObject* args)
{
  PyObject *v,*grad_w_dV,*v_X_grad_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &v,
                       &grad_w_dV,
                       &v_X_grad_w_dV))
    return NULL;
  calculateShape_X_weightedGradShape(SHAPE(v_X_grad_w_dV)[0],
                                     SHAPE(v_X_grad_w_dV)[1],
                                     SHAPE(v_X_grad_w_dV)[2],
                                     SHAPE(v_X_grad_w_dV)[3],
                                     SHAPE(v_X_grad_w_dV)[4],
                                     DDATA(v),
                                     DDATA(grad_w_dV),
                                     DDATA(v_X_grad_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateGradShape_X_weightedShape(PyObject* self, 
                                                PyObject* args)
{
  PyObject *grad_v,*w_dV,*grad_v_X_w_dV;
  if(!PyArg_ParseTuple(args,"OOO",
                       &grad_v,
                       &w_dV,
                       &grad_v_X_w_dV))
    return NULL;
  calculateGradShape_X_weightedShape(SHAPE(grad_v_X_w_dV)[0],
                                     SHAPE(grad_v_X_w_dV)[1],
                                     SHAPE(grad_v_X_w_dV)[2],
                                     SHAPE(grad_v_X_w_dV)[3],
                                     SHAPE(grad_v_X_w_dV)[4],
                                     DDATA(grad_v),
                                     DDATA(w_dV),
                                     DDATA(grad_v_X_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateGradShape_X_weightedGradShape(PyObject* self, 
                                                    PyObject* args)
{
  PyObject *grad_v,*grad_w_dV,*grad_v_X_grad_w_dV;  
  if(!PyArg_ParseTuple(args,"OOO",
                       &grad_v,
                       &grad_w_dV,
                       &grad_v_X_grad_w_dV))
    return NULL;
  calculateGradShape_X_weightedGradShape(SHAPE(grad_v_X_grad_w_dV)[0],
                                         SHAPE(grad_v_X_grad_w_dV)[1],
                                         SHAPE(grad_v_X_grad_w_dV)[2],
                                         SHAPE(grad_v_X_grad_w_dV)[3],
                                         SHAPE(grad_v_X_grad_w_dV)[4],
                                         DDATA(grad_v),
                                         DDATA(grad_w_dV),
                                         DDATA(grad_v_X_grad_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedShapeTrace(PyObject* self,
                                         PyObject* args)
{
  PyObject *dSR,*sqrt_det_g,*w,*w_dS;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dSR,
                       &sqrt_det_g,
                       &w,
                       &w_dS))
    return NULL;
  calculateWeightedShapeTrace(SHAPE(w_dS)[0],
                              SHAPE(w_dS)[1],
                              SHAPE(w_dS)[2],
                              SHAPE(w_dS)[3],
                              DDATA(dSR),
                              DDATA(sqrt_det_g),
                              DDATA(w),
                              DDATA(w_dS));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedPiolaShapeTrace(PyObject* self,
					      PyObject* args)
{
  PyObject *dSR,*sqrt_det_g,*w,*w_dS;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &dSR,
                       &sqrt_det_g,
                       &w,
                       &w_dS))
    return NULL;
  calculateWeightedPiolaShapeTrace(SHAPE(w_dS)[0],
				   SHAPE(w_dS)[1],
				   SHAPE(w_dS)[2],
				   SHAPE(w_dS)[3],
				   DDATA(dSR),
				   DDATA(sqrt_det_g),
				   DDATA(w),
				   DDATA(w_dS));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateShape_X_weightedShapeTrace(PyObject* self, 
                                                 PyObject* args)
{
  PyObject *v,*w_dS,*v_X_w_dS;
  if(!PyArg_ParseTuple(args,"OOO",
                       &v,
                       &w_dS,
                       &v_X_w_dS))
    return NULL;
  calculateShape_X_weightedShapeTrace(SHAPE(v_X_w_dS)[0],
                                      SHAPE(v_X_w_dS)[1],
                                      SHAPE(v_X_w_dS)[2],
                                      SHAPE(v_X_w_dS)[3],
                                      SHAPE(v_X_w_dS)[4],
                                      DDATA(v),
                                      DDATA(w_dS),
                                      DDATA(v_X_w_dS));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateGradShape_X_weightedShapeTrace(PyObject* self, 
                                                     PyObject* args)
{
  PyObject *grad_v,*w_dS,*grad_v_X_w_dS;
  if(!PyArg_ParseTuple(args,"OOO",
                       &grad_v,
                       &w_dS,
                       &grad_v_X_w_dS))
    return NULL;
  calculateGradShape_X_weightedShapeTrace(SHAPE(grad_v_X_w_dS)[0],
                                          SHAPE(grad_v_X_w_dS)[1],
                                          SHAPE(grad_v_X_w_dS)[2],
                                          SHAPE(grad_v_X_w_dS)[3],
                                          SHAPE(grad_v_X_w_dS)[4],
                                          SHAPE(grad_v_X_w_dS)[5],
                                          DDATA(grad_v),
                                          DDATA(w_dS),
                                          DDATA(grad_v_X_w_dS));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateWeightedShapeGlobalExteriorTrace(PyObject* self,
						       PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *dSR,*sqrt_det_g,*w,*w_dS;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
		       &elementBoundaryLocalElementBoundariesArray,
                       &dSR,
                       &sqrt_det_g,
                       &w,
                       &w_dS))
    return NULL;
  calculateWeightedShapeGlobalExteriorTrace(SHAPE(w_dS)[1],
					    SHAPE(w_dS)[2],
					    SHAPE(exteriorElementBoundariesArray)[0],
					    IDATA(exteriorElementBoundariesArray),
					    IDATA(elementBoundaryElementsArray),
					    IDATA(elementBoundaryLocalElementBoundariesArray),
					    DDATA(dSR),
					    DDATA(sqrt_det_g),
					    DDATA(w),
					    DDATA(w_dS));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cfemIntegralsCalculateShape_X_weightedShapeGlobalExteriorTrace(PyObject* self, 
							       PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *v,*w_dS,*v_X_w_dS;
  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
		       &elementBoundaryLocalElementBoundariesArray,
                       &v,
                       &w_dS,
                       &v_X_w_dS))
    return NULL;
  calculateShape_X_weightedShapeGlobalExteriorTrace(SHAPE(v_X_w_dS)[1],
						    SHAPE(v_X_w_dS)[2],
						    SHAPE(v_X_w_dS)[3],
						    SHAPE(exteriorElementBoundariesArray)[0],
						    IDATA(exteriorElementBoundariesArray),
						    IDATA(elementBoundaryElementsArray),
						    IDATA(elementBoundaryLocalElementBoundariesArray),
						    DDATA(v),
						    DDATA(w_dS),
						    DDATA(v_X_w_dS));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateGradShape_X_weightedShapeGlobalExteriorTrace(PyObject* self, 
								   PyObject* args)
{
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *grad_v,*w_dS,*grad_v_X_w_dS;
  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
		       &elementBoundaryLocalElementBoundariesArray,
                       &grad_v,
                       &w_dS,
                       &grad_v_X_w_dS))
    return NULL;
  calculateGradShape_X_weightedShapeGlobalExteriorTrace(SHAPE(grad_v_X_w_dS)[1],
							SHAPE(grad_v_X_w_dS)[2],
							SHAPE(grad_v_X_w_dS)[3],
							SHAPE(grad_v_X_w_dS)[4],
							SHAPE(exteriorElementBoundariesArray)[0],
							IDATA(exteriorElementBoundariesArray),
							IDATA(elementBoundaryElementsArray),
							IDATA(elementBoundaryLocalElementBoundariesArray),
							DDATA(grad_v),
							DDATA(w_dS),
							DDATA(grad_v_X_w_dS));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionValues(PyObject* self, 
                                                  PyObject* args)
{
  int nComponents=1;
  PyObject *l2g,*dof,*v,*u;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &l2g,
                       &dof,
                       &v,
                       &u))
    return NULL;
  if (ND(u) == 3)
    {
      nComponents = SHAPE(u)[2];
    }
  else
    {
      nComponents = 1;
    }
  calculateFiniteElementFunctionValues(SHAPE(v)[0],
                                       SHAPE(v)[1],
                                       SHAPE(v)[2],
                                       nComponents,
                                       IDATA(l2g),
                                       DDATA(dof),
                                       DDATA(v),
                                       DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionGradientValues(PyObject* self, 
                                                          PyObject* args)
{
  int nComponents=1,nSpace=1;
  PyObject *l2g,*dof,*grad_v,*grad_u;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &l2g,
                       &dof,
                       &grad_v,
                       &grad_u))
    return NULL;
  if (ND(grad_u) == 4)
    {
      nComponents = SHAPE(grad_u)[2];
      nSpace = SHAPE(grad_u)[3];
    }
  else
    {
      nComponents = 1;
      nSpace = SHAPE(grad_u)[2];
    }
  calculateFiniteElementFunctionGradientValues(SHAPE(grad_v)[0],
                                               SHAPE(grad_v)[1],
                                               SHAPE(grad_v)[2],
                                               nComponents,
                                               nSpace,
                                               IDATA(l2g),
                                               DDATA(dof),
                                               DDATA(grad_v),
                                               DDATA(grad_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionHessianValues(PyObject* self, 
							 PyObject* args)
{
  int nComponents=1,nSpace=1;
  PyObject *l2g,*dof,*Hessian_v,*Hessian_u;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &l2g,
                       &dof,
                       &Hessian_v,
                       &Hessian_u))
    return NULL;
  if (ND(Hessian_u) == 5)
    {
      nComponents = SHAPE(Hessian_u)[2];
      nSpace = SHAPE(Hessian_u)[3];
    }
  else
    {
      nComponents = 1;
      nSpace = SHAPE(Hessian_u)[2];
    }
  calculateFiniteElementFunctionHessianValues(SHAPE(Hessian_v)[0],
					      SHAPE(Hessian_v)[1],
					      SHAPE(Hessian_v)[2],
					      nComponents,
					      nSpace,
					      IDATA(l2g),
					      DDATA(dof),
					      DDATA(Hessian_v),
					      DDATA(Hessian_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionGradientTensorValues(PyObject* self, 
                                                                PyObject* args)
{
  int nComponents=1,nSpace=1;
  PyObject *l2g,*dof,*grad_v_X_grad_w_dV,*grad_u_X_grad_w_dV;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &l2g,
                       &dof,
                       &grad_v_X_grad_w_dV,
                       &grad_u_X_grad_w_dV)) 
    return NULL;
  if (ND(grad_u_X_grad_w_dV) == 6)
    {
      nComponents = SHAPE(grad_u_X_grad_w_dV)[3];
      nSpace = SHAPE(grad_u_X_grad_w_dV)[4];
    }
  else
    {
      nComponents = 1;
      nSpace = SHAPE(grad_u_X_grad_w_dV)[3];
    }
  calculateFiniteElementFunctionGradientTensorValues(SHAPE(grad_v_X_grad_w_dV)[0],
                                                     SHAPE(grad_v_X_grad_w_dV)[1],
                                                     SHAPE(grad_v_X_grad_w_dV)[2],
                                                     SHAPE(grad_v_X_grad_w_dV)[3],
                                                     nComponents,
                                                     nSpace,
                                                     IDATA(l2g),
                                                     DDATA(dof),
                                                     DDATA(grad_v_X_grad_w_dV),
                                                     DDATA(grad_u_X_grad_w_dV));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionValuesTrace(PyObject* self, 
                                                       PyObject* args)
{
  int nComponents=1;
  PyObject *l2g,*dof,*v,*u;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &l2g,
                       &dof,
                       &v,
                       &u))
    return NULL;
  if(ND(u) == 4)
    nComponents = SHAPE(u)[3];
  else
    nComponents = 1;
  calculateFiniteElementFunctionValuesTrace(SHAPE(v)[0],
                                            SHAPE(v)[1],
                                            SHAPE(v)[2],
                                            SHAPE(v)[3],
                                            nComponents,
                                            IDATA(l2g),
                                            DDATA(dof),
                                            DDATA(v),
                                            DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}


static PyObject*
cfemIntegralsCalculateFiniteElementFunctionGradientValuesTrace(PyObject* self, 
                                                               PyObject* args)
{
  int nComponents, nSpace;
  PyObject *l2g,*dof,*grad_v,*grad_u;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &l2g,
                       &dof,
                       &grad_v,
                       &grad_u))
    return NULL;
  if(ND(grad_u) == 5)
    {
      nComponents = SHAPE(grad_u)[3];
      nSpace = SHAPE(grad_u)[4];
    }
  else
    {
      nComponents = 1;
      nSpace = SHAPE(grad_u)[3];
    }
  calculateFiniteElementFunctionGradientValuesTrace(SHAPE(grad_v)[0],
                                                    SHAPE(grad_v)[1],
                                                    SHAPE(grad_v)[2],
                                                    SHAPE(grad_v)[3],
                                                    nComponents,
                                                    nSpace,
                                                    IDATA(l2g),
                                                    DDATA(dof),
                                                    DDATA(grad_v),
                                                    DDATA(grad_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionValuesGlobalExteriorTrace(PyObject* self, 
								     PyObject* args)
{
  int nComponents=1;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *l2g,*dof,*v,*u;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
		       &elementBoundaryLocalElementBoundariesArray,
                       &l2g,
                       &dof,
                       &v,
                       &u))
    return NULL;
  if(ND(u) == 4)
    nComponents = SHAPE(u)[3];
  else
    nComponents = 1;
  calculateFiniteElementFunctionValuesGlobalExteriorTrace(SHAPE(v)[1],
							  SHAPE(v)[2],
							  nComponents,
							  SHAPE(exteriorElementBoundariesArray)[0],
							  IDATA(exteriorElementBoundariesArray),
							  IDATA(elementBoundaryElementsArray),
							  IDATA(elementBoundaryLocalElementBoundariesArray),
							  IDATA(l2g),
							  DDATA(dof),
							  DDATA(v),
							  DDATA(u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFiniteElementFunctionGradientValuesGlobalExteriorTrace(PyObject* self, 
									     PyObject* args)
{
  int nComponents, nSpace;
  PyObject *exteriorElementBoundariesArray,
    *elementBoundaryElementsArray,
    *elementBoundaryLocalElementBoundariesArray,
    *l2g,*dof,*grad_v,*grad_u;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &exteriorElementBoundariesArray,
		       &elementBoundaryElementsArray,
		       &elementBoundaryLocalElementBoundariesArray,
                       &l2g,
                       &dof,
                       &grad_v,
                       &grad_u))
    return NULL;
  if(ND(grad_u) == 4)
    {
      nComponents = SHAPE(grad_u)[2];
      nSpace = SHAPE(grad_u)[3];
    }
  else
    {
      nComponents = 1;
      nSpace = SHAPE(grad_u)[2];
    }
  calculateFiniteElementFunctionGradientValuesGlobalExteriorTrace(SHAPE(grad_v)[1],
								  SHAPE(grad_v)[2],
								  nComponents,
								  nSpace,
								  SHAPE(exteriorElementBoundariesArray)[0],
								  IDATA(exteriorElementBoundariesArray),
								  IDATA(elementBoundaryElementsArray),
								  IDATA(elementBoundaryLocalElementBoundariesArray),
								  IDATA(l2g),
								  DDATA(dof),
								  DDATA(grad_v),
								  DDATA(grad_u));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalResidualFromElementResidual(PyObject* self, 
                                                     PyObject* args)
{
  int offset_r,stride_r;
  PyObject *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,*elementResidual,*globalResidual;
  if(!PyArg_ParseTuple(args,"iiOOOOO",
                       &offset_r,
                       &stride_r,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &elementResidual,
                       &globalResidual))
    return NULL;
  updateGlobalResidualFromElementResidual(SHAPE(elementResidual)[0],
                                          SHAPE(elementResidual)[1],
                                          offset_r,
                                          stride_r,
                                          IDATA(nFreeDOF_element_r),
                                          IDATA(freeLocal_r),
                                          IDATA(freeGlobal_r),
                                          DDATA(elementResidual),
                                          DDATA(globalResidual));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromElementJacobian_CSR(PyObject* self, 
                                                         PyObject* args)
{
  PyObject *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_ru,
    *elementJacobian,
    *globalJacobianCSR;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_ru,
                       &elementJacobian,
                       &globalJacobianCSR))
    return NULL;
  updateGlobalJacobianFromElementJacobian_CSR(SHAPE(elementJacobian)[0],
                                              SHAPE(elementJacobian)[1],
                                              SHAPE(elementJacobian)[2],
                                              IDATA(nFreeDOF_element_r),
                                              IDATA(freeLocal_r),
                                              IDATA(nFreeDOF_element_u),
                                              IDATA(freeLocal_u),
                                              IDATA(csrRowIndeces_ru),
                                              IDATA(csrColumnOffsets_ru),
                                              DDATA(elementJacobian),
                                              CSRVAL(globalJacobianCSR));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromElementJacobian_eb_CSR(PyObject* self, 
                                                         PyObject* args)
{
  PyObject *elementNeighbors,
    *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_eb_ru,
    *elementJacobian_eb,
    *globalJacobianCSR;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                       &elementNeighbors,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_eb_ru,
                       &elementJacobian_eb,
                       &globalJacobianCSR))
    return NULL;
  updateGlobalJacobianFromElementJacobian_eb_CSR(IDATA(elementNeighbors),
                                                 SHAPE(elementJacobian_eb)[0],
                                                 SHAPE(elementJacobian_eb)[1],
                                                 SHAPE(elementJacobian_eb)[2],
                                                 SHAPE(elementJacobian_eb)[3],
                                                 IDATA(nFreeDOF_element_r),
                                                 IDATA(freeLocal_r),
                                                 IDATA(nFreeDOF_element_u),
                                                 IDATA(freeLocal_u),
                                                 IDATA(csrRowIndeces_ru),
                                                 IDATA(csrColumnOffsets_eb_ru),
                                                 DDATA(elementJacobian_eb),
                                                 CSRVAL(globalJacobianCSR));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromElementJacobian_dense(PyObject* self, 
                                                           PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  PyObject *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,*nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,*elementJacobian,*jacobian;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOO",
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
		       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
		       &freeGlobal_u,
                       &elementJacobian,
                       &jacobian))
    return NULL;
  updateGlobalJacobianFromElementJacobian_dense(SHAPE(elementJacobian)[0],
                                                SHAPE(elementJacobian)[1],
                                                SHAPE(elementJacobian)[2],
                                                offset_r,
                                                stride_r,
                                                offset_u,
                                                stride_u,
						nFreeVDOF_global,
						IDATA(nFreeDOF_element_r),
						IDATA(freeLocal_r),
						IDATA(freeGlobal_r),
						IDATA(nFreeDOF_element_u),
						IDATA(freeLocal_u),
						IDATA(freeGlobal_u),
						DDATA(elementJacobian),
						DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromElementJacobian_eb_dense(PyObject* self, 
                                                           PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  PyObject *elementNeighbors,*nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,*nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,*elementJacobian_eb,*jacobian;
  if(!PyArg_ParseTuple(args,"OiiiiiOOOOOOOO",
                       &elementNeighbors,
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
		       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
		       &freeGlobal_u,
                       &elementJacobian_eb,
                       &jacobian))
    return NULL;
  updateGlobalJacobianFromElementJacobian_eb_dense(IDATA(elementNeighbors),
                                                   SHAPE(elementJacobian_eb)[0],
                                                   SHAPE(elementJacobian_eb)[1],
                                                   SHAPE(elementJacobian_eb)[2],
                                                   SHAPE(elementJacobian_eb)[3],
                                                   offset_r,
                                                   stride_r,
                                                   offset_u,
                                                   stride_u,
                                                   nFreeVDOF_global,
                                                   IDATA(nFreeDOF_element_r),
                                                   IDATA(freeLocal_r),
                                                   IDATA(freeGlobal_r),
                                                   IDATA(nFreeDOF_element_u),
                                                   IDATA(freeLocal_u),
                                                   IDATA(freeGlobal_u),
                                                   DDATA(elementJacobian_eb),
                                                   DDATA(jacobian));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateFlowVelocity(PyObject* self, 
                                   PyObject* args)
{
  PyObject *f,*a,*grad_phi,*v;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &f,
                       &a,
                       &grad_phi,
                       &v))
    return NULL;
  calculateFlowVelocity(SHAPE(f)[0],
                        SHAPE(f)[1],
                        SHAPE(f)[2],
                        DDATA(f),
                        DDATA(a),
                        DDATA(grad_phi),
                        DDATA(v));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateAddJacobian_CSR(PyObject* self, 
                                   PyObject* args)
{
  int jacIndex;
  double  val;
  PyObject *jac;
  
  if(!PyArg_ParseTuple(args,"idO",
                       &jacIndex,
                       &val,
                       &jac))
    return NULL;
  updateAddJacobian_CSR(jacIndex,val,CSRVAL(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsZeroJacobian_CSR(PyObject* self, 
                              PyObject* args)
{
  int nNonzeros;
  PyObject *jac;
  if(!PyArg_ParseTuple(args,"iO",
                       &nNonzeros,
                       &jac))
    return NULL;
  zeroJacobian_CSR(nNonzeros,CSRVAL(jac));
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(PyObject* self, 
                                                                             PyObject* args)
{
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_eb_ru,
    *elementBoundaryFluxJacobian,
    *w_dS,
    *jac;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_eb_ru,
                       &elementBoundaryFluxJacobian,
                       &w_dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR(SHAPE(interiorElementBoundaries)[0],
                                                                  SHAPE(w_dS)[1],
                                                                  SHAPE(w_dS)[2],
                                                                  SHAPE(w_dS)[3],
                                                                  SHAPE(elementBoundaryFluxJacobian)[3],
                                                                  IDATA(interiorElementBoundaries),
                                                                  IDATA(elementBoundaryElements),
                                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                                  IDATA(nFreeDOF_element_r),
                                                                  IDATA(freeLocal_r),
                                                                  IDATA(nFreeDOF_element_u),
                                                                  IDATA(freeLocal_u),
                                                                  IDATA(csrRowIndeces_ru),
                                                                  IDATA(csrColumnOffsets_eb_ru),
                                                                  DDATA(elementBoundaryFluxJacobian),
                                                                  DDATA(w_dS),
                                                                  CSRVAL(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_CSR(PyObject* self, 
                                                                             PyObject* args)
{
  PyObject *elementNeighbors,
    *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_eb_eNebN_ru,
    *elementBoundaryFluxJacobian_eb,
    *w_dS,
    *jac;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
                       &elementNeighbors,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_eb_eNebN_ru,
                       &elementBoundaryFluxJacobian_eb,
                       &w_dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_CSR(IDATA(elementNeighbors),
                                                                     SHAPE(interiorElementBoundaries)[0],
                                                                  SHAPE(w_dS)[1],
                                                                  SHAPE(w_dS)[2],
                                                                  SHAPE(w_dS)[3],
                                                                  SHAPE(elementBoundaryFluxJacobian_eb)[4],
                                                                  IDATA(interiorElementBoundaries),
                                                                  IDATA(elementBoundaryElements),
                                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                                  IDATA(nFreeDOF_element_r),
                                                                  IDATA(freeLocal_r),
                                                                  IDATA(nFreeDOF_element_u),
                                                                  IDATA(freeLocal_u),
                                                                  IDATA(csrRowIndeces_ru),
                                                                  IDATA(csrColumnOffsets_eb_eNebN_ru),
                                                                  DDATA(elementBoundaryFluxJacobian_eb),
                                                                  DDATA(w_dS),
                                                                  CSRVAL(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_CSR(PyObject* self, 
										    PyObject* args)
{
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_eb_ru,
    *elementBoundaryFluxJacobian_2sided,
    *w_dS,
    *jac;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_eb_ru,
                       &elementBoundaryFluxJacobian_2sided,
                       &w_dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_CSR(SHAPE(interiorElementBoundaries)[0],
									 SHAPE(w_dS)[1],
									 SHAPE(w_dS)[2],
									 SHAPE(w_dS)[3],
									 SHAPE(elementBoundaryFluxJacobian_2sided)[4],
									 IDATA(interiorElementBoundaries),
									 IDATA(elementBoundaryElements),
									 IDATA(elementBoundaryLocalElementBoundaries),
									 IDATA(nFreeDOF_element_r),
									 IDATA(freeLocal_r),
									 IDATA(nFreeDOF_element_u),
									 IDATA(freeLocal_u),
									 IDATA(csrRowIndeces_ru),
									 IDATA(csrColumnOffsets_eb_ru),
									 DDATA(elementBoundaryFluxJacobian_2sided),
									 DDATA(w_dS),
									 CSRVAL(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(PyObject* self, 
                                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_eb_ru,
    *elementBoundaryFluxJacobian,
    *w_dS,
    *jac;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_eb_ru,
                       &elementBoundaryFluxJacobian,
                       &w_dS,
                       &jac))
    return NULL;
  if (ND(w_dS) > 3)
    {
      assert(ND(elementBoundaryFluxJacobian) == 4);
      updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR(SHAPE(exteriorElementBoundaries)[0],
								      SHAPE(w_dS)[1],
								      SHAPE(w_dS)[2],
								      SHAPE(w_dS)[3],
								      SHAPE(elementBoundaryFluxJacobian)[3],
								      IDATA(exteriorElementBoundaries),
								      IDATA(elementBoundaryElements),
								      IDATA(elementBoundaryLocalElementBoundaries),
								      IDATA(nFreeDOF_element_r),
								      IDATA(freeLocal_r),
								      IDATA(nFreeDOF_element_u),
								      IDATA(freeLocal_u),
								      IDATA(csrRowIndeces_ru),
								      IDATA(csrColumnOffsets_eb_ru),
								      DDATA(elementBoundaryFluxJacobian),
								      DDATA(w_dS),
								      CSRVAL(jac));
    }
  else
    {
      assert(ND(elementBoundaryFluxJacobian) == 3);
      updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_CSR(SHAPE(exteriorElementBoundaries)[0],
									    SHAPE(w_dS)[1],
									    SHAPE(w_dS)[2],
									    SHAPE(elementBoundaryFluxJacobian)[2],
									    IDATA(exteriorElementBoundaries),
									    IDATA(elementBoundaryElements),
									    IDATA(elementBoundaryLocalElementBoundaries),
									    IDATA(nFreeDOF_element_r),
									    IDATA(freeLocal_r),
									    IDATA(nFreeDOF_element_u),
									    IDATA(freeLocal_u),
									    IDATA(csrRowIndeces_ru),
									    IDATA(csrColumnOffsets_eb_ru),
									    DDATA(elementBoundaryFluxJacobian),
									    DDATA(w_dS),
									    CSRVAL(jac));
      
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_CSR(PyObject* self, 
                                                                             PyObject* args)
{
  PyObject *elementNeighbors,*exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,
    *nFreeDOF_element_u,*freeLocal_u,
    *csrRowIndeces_ru,*csrColumnOffsets_eb_eNebN_ru,
    *elementBoundaryFluxJacobian_eb,
    *w_dS,
    *jac;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOO",
                       &elementNeighbors,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &csrRowIndeces_ru,
                       &csrColumnOffsets_eb_eNebN_ru,
                       &elementBoundaryFluxJacobian_eb,
                       &w_dS,
                       &jac))
    return NULL;
  if (ND(w_dS) > 3)
    {
      updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_CSR(IDATA(elementNeighbors),
									 SHAPE(exteriorElementBoundaries)[0],
									 SHAPE(w_dS)[1],
									 SHAPE(w_dS)[2],
									 SHAPE(w_dS)[3],
									 SHAPE(elementBoundaryFluxJacobian_eb)[4],
									 IDATA(exteriorElementBoundaries),
									 IDATA(elementBoundaryElements),
									 IDATA(elementBoundaryLocalElementBoundaries),
									 IDATA(nFreeDOF_element_r),
									 IDATA(freeLocal_r),
									 IDATA(nFreeDOF_element_u),
									 IDATA(freeLocal_u),
									 IDATA(csrRowIndeces_ru),
									 IDATA(csrColumnOffsets_eb_eNebN_ru),
									 DDATA(elementBoundaryFluxJacobian_eb),
									 DDATA(w_dS),
									 CSRVAL(jac));
    }
  else
    {
      updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_eb_CSR(IDATA(elementNeighbors),
									       SHAPE(exteriorElementBoundaries)[0],
									       SHAPE(elementBoundaryFluxJacobian_eb)[2],
									       SHAPE(w_dS)[1],
									       SHAPE(w_dS)[2],
									       SHAPE(elementBoundaryFluxJacobian_eb)[4],
									       IDATA(exteriorElementBoundaries),
									       IDATA(elementBoundaryElements),
									       IDATA(elementBoundaryLocalElementBoundaries),
									       IDATA(nFreeDOF_element_r),
									       IDATA(freeLocal_r),
									       IDATA(nFreeDOF_element_u),
									       IDATA(freeLocal_u),
									       IDATA(csrRowIndeces_ru),
									       IDATA(csrColumnOffsets_eb_eNebN_ru),
									       DDATA(elementBoundaryFluxJacobian_eb),
									       DDATA(w_dS),
									       CSRVAL(jac));
     }
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(PyObject* self, 
                                                                               PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *elementBoundaryFluxJacobian,*w_dS,*jac;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOOOO",
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &elementBoundaryFluxJacobian,
                       &w_dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense(SHAPE(interiorElementBoundaries)[0],
                                                                    SHAPE(w_dS)[1],
								    SHAPE(w_dS)[2],
                                                                    SHAPE(w_dS)[3],
                                                                    SHAPE(elementBoundaryFluxJacobian)[3],
                                                                    offset_r,
                                                                    stride_r,
                                                                    offset_u,
                                                                    stride_u,
								    nFreeVDOF_global,
                                                                    IDATA(interiorElementBoundaries),
                                                                    IDATA(elementBoundaryElements),
                                                                    IDATA(elementBoundaryLocalElementBoundaries),
                                                                    IDATA(nFreeDOF_element_r),
                                                                    IDATA(nFreeDOF_element_u),
                                                                    IDATA(freeLocal_r),
                                                                    IDATA(freeGlobal_r),
                                                                    IDATA(freeLocal_u),
                                                                    IDATA(freeGlobal_u),
                                                                    DDATA(elementBoundaryFluxJacobian),
                                                                    DDATA(w_dS),
                                                                    DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_dense(PyObject* self, 
                                                                               PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global,nElements_global;
  PyObject *elementNeighbors,*interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *elementBoundaryFluxJacobian_eb,*w_dS,*jac;
  if(!PyArg_ParseTuple(args,"OiiiiiiOOOOOOOOOOOO",
                       &elementNeighbors,
                       &nElements_global,
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &elementBoundaryFluxJacobian_eb,
                       &w_dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_dense(IDATA(elementNeighbors),
                                                                       nElements_global,
                                                                       SHAPE(interiorElementBoundaries)[0],
                                                                       SHAPE(w_dS)[1],
                                                                       SHAPE(w_dS)[2],
                                                                       SHAPE(w_dS)[3],
                                                                       SHAPE(elementBoundaryFluxJacobian_eb)[4],
                                                                       offset_r,
                                                                       stride_r,
                                                                       offset_u,
                                                                       stride_u,
                                                                       nFreeVDOF_global,
                                                                       IDATA(interiorElementBoundaries),
                                                                       IDATA(elementBoundaryElements),
                                                                       IDATA(elementBoundaryLocalElementBoundaries),
                                                                       IDATA(nFreeDOF_element_r),
                                                                       IDATA(nFreeDOF_element_u),
                                                                       IDATA(freeLocal_r),
                                                                       IDATA(freeGlobal_r),
                                                                       IDATA(freeLocal_u),
                                                                       IDATA(freeGlobal_u),
                                                                       DDATA(elementBoundaryFluxJacobian_eb),
                                                                       DDATA(w_dS),
                                                                       DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_dense(PyObject* self, 
										      PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *elementBoundaryFluxJacobian_2sided,*w_dS,*jac;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOOOO",
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &elementBoundaryFluxJacobian_2sided,
                       &w_dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_dense(SHAPE(interiorElementBoundaries)[0],
									   SHAPE(w_dS)[1],
									   SHAPE(w_dS)[2],
									   SHAPE(w_dS)[3],
									   SHAPE(elementBoundaryFluxJacobian_2sided)[4],
									   offset_r,
									   stride_r,
									   offset_u,
									   stride_u,
									   nFreeVDOF_global,
									   IDATA(interiorElementBoundaries),
									   IDATA(elementBoundaryElements),
									   IDATA(elementBoundaryLocalElementBoundaries),
									   IDATA(nFreeDOF_element_r),
									   IDATA(nFreeDOF_element_u),
									   IDATA(freeLocal_r),
									   IDATA(freeGlobal_r),
									   IDATA(freeLocal_u),
									   IDATA(freeGlobal_u),
									   DDATA(elementBoundaryFluxJacobian_2sided),
									   DDATA(w_dS),
									   DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(PyObject* self, 
                                                                               PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *elementBoundaryFluxJacobian,*w_dS,*jac;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOOOO",
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &elementBoundaryFluxJacobian,
                       &w_dS,
                       &jac))
    return NULL;
  if (ND(w_dS) > 3)
    {
      assert(ND(elementBoundaryFluxJacobian) == 4);
      updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense(SHAPE(exteriorElementBoundaries)[0],
									SHAPE(w_dS)[1],
									SHAPE(w_dS)[2],
									SHAPE(w_dS)[3],
									SHAPE(elementBoundaryFluxJacobian)[3],
									offset_r,
									stride_r,
									offset_u,
									stride_u,
									nFreeVDOF_global,
									IDATA(exteriorElementBoundaries),
									IDATA(elementBoundaryElements),
									IDATA(elementBoundaryLocalElementBoundaries),
									IDATA(nFreeDOF_element_r),
									IDATA(nFreeDOF_element_u),
									IDATA(freeLocal_r),
									IDATA(freeGlobal_r),
									IDATA(freeLocal_u),
									IDATA(freeGlobal_u),
									DDATA(elementBoundaryFluxJacobian),
									DDATA(w_dS),
									DDATA(jac));
    }
  else
    {
      assert(ND(elementBoundaryFluxJacobian) == 3);
      updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_dense(SHAPE(exteriorElementBoundaries)[0],
									      SHAPE(w_dS)[1],
									      SHAPE(w_dS)[2],
									      SHAPE(elementBoundaryFluxJacobian)[2],
									      offset_r,
									      stride_r,
									      offset_u,
									      stride_u,
									      nFreeVDOF_global,
									      IDATA(exteriorElementBoundaries),
									      IDATA(elementBoundaryElements),
									      IDATA(elementBoundaryLocalElementBoundaries),
									      IDATA(nFreeDOF_element_r),
									      IDATA(nFreeDOF_element_u),
									      IDATA(freeLocal_r),
									      IDATA(freeGlobal_r),
									      IDATA(freeLocal_u),
									      IDATA(freeGlobal_u),
									      DDATA(elementBoundaryFluxJacobian),
									      DDATA(w_dS),
									      DDATA(jac));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_dense(PyObject* self, 
                                                                               PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global,nElements_global;
  PyObject *elementNeighbors,*exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *elementBoundaryFluxJacobian_eb,*w_dS,*jac;
  if(!PyArg_ParseTuple(args,"OiiiiiiOOOOOOOOOOOO",
                       &elementNeighbors,
                       &nElements_global,
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &elementBoundaryFluxJacobian_eb,
                       &w_dS,
                       &jac))
    return NULL;
  if (ND(w_dS) > 3)
    {
      updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_dense(IDATA(elementNeighbors),
									   nElements_global,
									   SHAPE(exteriorElementBoundaries)[0],
									   SHAPE(w_dS)[1],
									   SHAPE(w_dS)[2],
									   SHAPE(w_dS)[3],
									   SHAPE(elementBoundaryFluxJacobian_eb)[4],
									   offset_r,
									   stride_r,
									   offset_u,
									   stride_u,
									   nFreeVDOF_global,
									   IDATA(exteriorElementBoundaries),
									   IDATA(elementBoundaryElements),
									   IDATA(elementBoundaryLocalElementBoundaries),
									   IDATA(nFreeDOF_element_r),
									   IDATA(nFreeDOF_element_u),
									   IDATA(freeLocal_r),
									   IDATA(freeGlobal_r),
									   IDATA(freeLocal_u),
									   IDATA(freeGlobal_u),
									   DDATA(elementBoundaryFluxJacobian_eb),
									   DDATA(w_dS),
									   DDATA(jac));
    }
  else
    {
      updateGlobalJacobianFromGlobalExteriorElementBoundaryFluxJacobian_eb_dense(IDATA(elementNeighbors),
										 nElements_global,
										 SHAPE(exteriorElementBoundaries)[0],
										 SHAPE(elementBoundaryFluxJacobian_eb)[2],
										 SHAPE(w_dS)[1],
										 SHAPE(w_dS)[2],
										 SHAPE(elementBoundaryFluxJacobian_eb)[4],
										 offset_r,
										 stride_r,
										 offset_u,
										 stride_u,
										 nFreeVDOF_global,
										 IDATA(exteriorElementBoundaries),
										 IDATA(elementBoundaryElements),
										 IDATA(elementBoundaryLocalElementBoundaries),
										 IDATA(nFreeDOF_element_r),
										 IDATA(nFreeDOF_element_u),
										 IDATA(freeLocal_r),
										 IDATA(freeGlobal_r),
										 IDATA(freeLocal_u),
										 IDATA(freeGlobal_u),
										 DDATA(elementBoundaryFluxJacobian_eb),
										 DDATA(w_dS),
										 DDATA(jac));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryFlux(PyObject* self, 
                                               PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *flux,
    *w,
    *residual;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &flux,
                       &w,
                       &residual))
    return NULL;
  if (ND(w) >= 4)
    updateExteriorElementBoundaryFlux(SHAPE(exteriorElementBoundaries)[0],
				      SHAPE(w)[1],
				      SHAPE(w)[2],
				      SHAPE(w)[3],
				      IDATA(exteriorElementBoundaries),
				      IDATA(elementBoundaryElements),
				      IDATA(elementBoundaryLocalElementBoundaries),
				      DDATA(flux),
				      DDATA(w),
				      DDATA(residual));
  else
    updateGlobalExteriorElementBoundaryFlux(SHAPE(exteriorElementBoundaries)[0],
					    SHAPE(w)[1],
					    SHAPE(w)[2],
					    IDATA(exteriorElementBoundaries),
					    IDATA(elementBoundaryElements),
					    IDATA(elementBoundaryLocalElementBoundaries),
					    DDATA(flux),
					    DDATA(w),
					    DDATA(residual));
 
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryStressFlux(PyObject* self, 
						     PyObject* args)
{
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *stressFlux,
    *w,
    *residual;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &stressFlux,
                       &w,
                       &residual))
    return NULL;
  updateGlobalExteriorElementBoundaryStressFlux(SHAPE(exteriorElementBoundaries)[0],
						SHAPE(w)[1],
						SHAPE(w)[2],
						IDATA(exteriorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						DDATA(stressFlux),
						DDATA(w),
						DDATA(residual));
  
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryFlux(PyObject* self, 
                                               PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*flux,*w,*residual;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &flux,
                       &w,
                       &residual))
    return NULL;
  updateInteriorElementBoundaryFlux(SHAPE(interiorElementBoundaries)[0],
                                    SHAPE(w)[1],
                                    SHAPE(w)[2],
                                    SHAPE(w)[3],
                                    IDATA(interiorElementBoundaries),
                                    IDATA(elementBoundaryElements),
                                    IDATA(elementBoundaryLocalElementBoundaries),
                                    DDATA(flux),
                                    DDATA(w),
                                    DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}




static PyObject*
cfemIntegralsUpdateTwoSidedInteriorElementBoundaryFlux(PyObject* self, 
						       PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*flux,*w,*residual;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &flux,
                       &w,
                       &residual))
    return NULL;
  updateInteriorTwoSidedElementBoundaryFlux(SHAPE(interiorElementBoundaries)[0],
					    SHAPE(w)[1],
					    SHAPE(w)[2],
					    SHAPE(w)[3],
					    IDATA(interiorElementBoundaries),
					    IDATA(elementBoundaryElements),
					    IDATA(elementBoundaryLocalElementBoundaries),
					    DDATA(flux),
					    DDATA(w),
					    DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}



static PyObject*
cfemIntegralsCalculateInteriorElementBoundaryVelocities(PyObject* self, 
                                                        PyObject* args)
{
  /*   int nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace; */
  /*   PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*m,*a,*grad_phi,*f,*vAverage,*vJump,*mAverage,*mJump; */
  /*   if(!PyArg_ParseTuple(args,"iiiiOOOOOOOOOOO", */
  /*                        &nQuadraturePoints_elementBoundary, */
  /*                        &nElementBoundaries_element, */
  /*                        &nSpace, */
  /*                        &interiorElementBoundaries, */
  /*                        &elementBoundaryElements, */
  /*                        &elementBoundaryLocalElementBoundaries, */
  /*                        &m, */
  /*                        &a, */
  /* 		       &grad_phi, */
  /* 		       &f, */
  /* 		       &vAverage, */
  /* 		       &vJump, */
  /* 		       &mAverage, */
  /* 		       &mJump)) */
  /*     return NULL; */
  /*   calculateInteriorElementBoundaryVelocities(SHAPE(interiorElementBoundaries)[0], */
  /* 					     nElementBoundaries_element, */
  /* 					     nQuadraturePoints_elementBoundary, */
  /* 					     nSpace, */
  /* 					     IDATA(interiorElementBoundaries), */
  /* 					     IDATA(elementBoundaryElements), */
  /* 					     IDATA(elementBoundaryLocalElementBoundaries), */
  /* 					     DDATA(m), */
  /* 					     DDATA(a), */
  /* 					     DDATA(grad_phi), */
  /* 					     DDATA(f), */
  /* 					     DDATA(vAverage), */
  /* 					     DDATA(vJump), */
  /* 					     DDATA(mAverage), */
  /* 					     DDATA(mJump)); */
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateExteriorElementBoundaryVelocities(PyObject* self, 
                                                        PyObject* args)
{
  /*   int nQuadraturePoints_elementBoundary,nElementBoundaries_element,nSpace; */
  /*   PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*m,*a,*grad_phi,*f,*vAverage,*vJump,*mAverage,*mJump; */
  /*   if(!PyArg_ParseTuple(args,"iiiOOOOOOOOOOO", */
  /*                        &nQuadraturePoints_elementBoundary, */
  /*                        &nElementBoundaries_element, */
  /*                        &nSpace, */
  /*                        &exteriorElementBoundaries, */
  /*                        &elementBoundaryElements, */
  /*                        &elementBoundaryLocalElementBoundaries, */
  /*                        &m, */
  /*                        &a, */
  /* 		       &grad_phi, */
  /* 		       &f, */
  /* 		       &vAverage, */
  /* 		       &vJump, */
  /* 		       &mAverage, */
  /* 		       &mJump)) */
  /*     return NULL; */
  /*   calculateExteriorElementBoundaryVelocities(SHAPE(exteriorElementBoundaries)[0], */
  /* 					     nElementBoundaries_element, */
  /* 					     nQuadraturePoints_elementBoundary, */
  /* 					     nSpace, */
  /* 					     IDATA(exteriorElementBoundaries), */
  /* 					     IDATA(elementBoundaryElements), */
  /* 					     IDATA(elementBoundaryLocalElementBoundaries), */
  /* 					     DDATA(m), */
  /* 					     DDATA(a), */
  /* 					     DDATA(grad_phi), */
  /* 					     DDATA(f), */
  /* 					     DDATA(vAverage), */
  /* 					     DDATA(vJump), */
  /* 					     DDATA(mAverage), */
  /* 					     DDATA(mJump)); */
  Py_INCREF(Py_None);
  return Py_None;
}




static PyObject*
cfemIntegralsWriteDOF(PyObject* self, 
                      PyObject* args)
{
  int j,nDOF,nComponents_DOF,component;
  const  char* format;
  double* dofptr;
  FILE* fileptr;
  PyObject *dof,*file;
  if(!PyArg_ParseTuple(args,"iiisOO",
                       &nDOF,
                       &nComponents_DOF,
                       &component,
                       &format,
                       &dof,
                       &file))
    return NULL;
  dofptr = DDATA(dof);
  fileptr = PyFile_AsFile(file);
  for(j=0;j<nDOF;j++)
    fprintf(fileptr,format,dofptr[j*nComponents_DOF+component]);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsWriteDOF_ZEROS(PyObject* self, 
                      PyObject* args)
{
  int j,nDOF,nComponents_DOF,component;
  const  char* format;
  FILE* fileptr;
  PyObject *file;
  if(!PyArg_ParseTuple(args,"iiisO",
                       &nDOF,
                       &nComponents_DOF,
                       &component,
                       &format,
                       &file))
    return NULL;
  fileptr = PyFile_AsFile(file);
  for(j=0;j<nDOF;j++)
    fprintf(fileptr,format,0.0);
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateDimensionlessNumbersADR(PyObject* self,
					      PyObject* args)
{
  PyObject *elementDiameter,*df,*a,*dphi,*dr,*dmt,
    *pe,*cfl;
  int nElements_global,nQuadraturePoints_element,nSpace_global;
  int computeDiffusiveTimeStepLimit=0;
  if(!PyArg_ParseTuple(args,"iiiOOOOOOOO",
		       &nElements_global,
		       &nQuadraturePoints_element,
		       &nSpace_global,
                       &elementDiameter,
                       &df,
                       &a,
                       &dphi,
                       &dr,
                       &dmt,
                       &pe,
                       &cfl))
    return NULL;
  calculateDimensionlessNumbersADR(nElements_global,
				   nQuadraturePoints_element,
				   nSpace_global,
				   computeDiffusiveTimeStepLimit,
                                   DDATA(elementDiameter),
                                   DDATA(df),
                                   DDATA(a),
				   DDATA(dphi),
                                   DDATA(dr),
                                   DDATA(dmt),
                                   DDATA(pe),
                                   DDATA(cfl));
                                   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateDimensionlessNumbersADR_sd(PyObject* self,
                                                 PyObject* args)
{
  PyObject *elementDiameter,*df,*a,*dphi,*dr,*dmt,
    *pe,*cfl,*rowptr,*colind;
  int nElements_global,nQuadraturePoints_element,nSpace_global;
  int computeDiffusiveTimeStepLimit=0;
  if(!PyArg_ParseTuple(args,"iiiOOOOOOOOOO",
		       &nElements_global,
		       &nQuadraturePoints_element,
		       &nSpace_global,
                       &rowptr,
                       &colind,
                       &elementDiameter,
                       &df,
                       &a,
                       &dphi,
                       &dr,
                       &dmt,
                       &pe,
                       &cfl))
    return NULL;
  calculateDimensionlessNumbersADR_sd(nElements_global,
                                      nQuadraturePoints_element,
                                      nSpace_global,
				      computeDiffusiveTimeStepLimit,
                                      IDATA(rowptr),
                                      IDATA(colind),
                                      DDATA(elementDiameter),
                                      DDATA(df),
                                      DDATA(a),
                                      DDATA(dphi),
                                      DDATA(dr),
                                      DDATA(dmt),
                                      DDATA(pe),
                                      DDATA(cfl));
                                   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateCFLADR(PyObject* self,
			     PyObject* args)
{
  PyObject *elementDiameter,*dm,*df,*cfl;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &elementDiameter,
		       &dm,
                       &df,
                       &cfl))
    return NULL;
  calculateCFLADR(SHAPE(df)[0],
		  SHAPE(df)[1],
		  SHAPE(df)[2],
		  DDATA(elementDiameter),
		  DDATA(dm),
		  DDATA(df),
		  DDATA(cfl));
                                   
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cfemIntegralsCalculateCFLADR2speeds(PyObject* self,
				    PyObject* args)
{
  PyObject *elementDiameter,*dm,*df1,*df2,*cfl;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &elementDiameter,
		       &dm,
                       &df1,
		       &df2,
                       &cfl))
    return NULL;
  calculateCFLADR2speeds(SHAPE(df1)[0],
			 SHAPE(df1)[1],
			 SHAPE(df1)[2],
			 DDATA(elementDiameter),
			 DDATA(dm),
			 DDATA(df1),
			 DDATA(df2),
			 DDATA(cfl));
                                   
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateIntegrationWeights(PyObject* self,
					 PyObject* args)
{
  
  PyObject *abs_det_J,*referenceWeights,*physicalWeights;
  if(!PyArg_ParseTuple(args,"OOO",
		       &abs_det_J,
                       &referenceWeights,
		       &physicalWeights))    
    return NULL;

  calculateIntegrationWeights(SHAPE(abs_det_J)[0],
			      SHAPE(abs_det_J)[1],
			      DDATA(abs_det_J),
			      DDATA(referenceWeights),
			      DDATA(physicalWeights));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateElementBoundaryIntegrationWeights(PyObject* self,
							PyObject* args)
{
  
  PyObject *sqrt_det_g,*referenceWeights,*physicalWeights;
  if(!PyArg_ParseTuple(args,"OOO",
		       &sqrt_det_g,
                       &referenceWeights,
		       &physicalWeights))    
    return NULL;

  calculateElementBoundaryIntegrationWeights(SHAPE(sqrt_det_g)[0],
					     SHAPE(sqrt_det_g)[1],
					     SHAPE(sqrt_det_g)[2],
					     DDATA(sqrt_det_g),
					     DDATA(referenceWeights),
					     DDATA(physicalWeights));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryAdvectiveVelocity(PyObject* self, 
                                                            PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*f,*velocity;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &f,
                       &velocity))
    return NULL;
  updateInteriorElementBoundaryAdvectiveVelocity(SHAPE(interiorElementBoundaries)[0],
                                                 SHAPE(f)[1],
                                                 SHAPE(f)[2],
                                                 SHAPE(f)[3],
                                                 IDATA(interiorElementBoundaries),
                                                 IDATA(elementBoundaryElements),
                                                 IDATA(elementBoundaryLocalElementBoundaries),
                                                 DDATA(f),
                                                 DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryAdvectiveVelocity(PyObject* self, 
                                                            PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*f,*velocity;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &f,
                       &velocity))
    return NULL;
  if (ND(f) > 3)
    {
      assert(ND(velocity) == ND(f));
      updateExteriorElementBoundaryAdvectiveVelocity(SHAPE(exteriorElementBoundaries)[0],
						     SHAPE(f)[1],
						     SHAPE(f)[2],
						     SHAPE(f)[3],
						     IDATA(exteriorElementBoundaries),
						     IDATA(elementBoundaryElements),
						     IDATA(elementBoundaryLocalElementBoundaries),
						     DDATA(f),
						     DDATA(velocity));
    }
  else
    {
      assert(ND(velocity) == ND(f));
      updateGlobalExteriorElementBoundaryAdvectiveVelocity(SHAPE(exteriorElementBoundaries)[0],
							   SHAPE(f)[1],
							   SHAPE(f)[2],
							   IDATA(exteriorElementBoundaries),
							   IDATA(elementBoundaryElements),
							   IDATA(elementBoundaryLocalElementBoundaries),
							   DDATA(f),
							   DDATA(velocity));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryDiffusiveVelocity(PyObject* self, 
                                                            PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*a,*grad_phi,*velocity;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &a,
                       &grad_phi,
                       &velocity))
    return NULL;
  updateInteriorElementBoundaryDiffusiveVelocity(SHAPE(interiorElementBoundaries)[0],
                                                 SHAPE(grad_phi)[1],
                                                 SHAPE(grad_phi)[2],
                                                 SHAPE(grad_phi)[3],
                                                 IDATA(interiorElementBoundaries),
                                                 IDATA(elementBoundaryElements),
                                                 IDATA(elementBoundaryLocalElementBoundaries),
                                                 DDATA(a),
                                                 DDATA(grad_phi),
                                                 DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryDiffusiveVelocity_sd(PyObject* self, 
                                                            PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*a,*grad_phi,*velocity,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
                       &rowptr,
                       &colind,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &a,
                       &grad_phi,
                       &velocity))
    return NULL;
  updateInteriorElementBoundaryDiffusiveVelocity_sd(SHAPE(interiorElementBoundaries)[0],
                                                    SHAPE(grad_phi)[1],
                                                    SHAPE(grad_phi)[2],
                                                    SHAPE(grad_phi)[3],
                                                    IDATA(rowptr),
                                                    IDATA(colind),
                                                    IDATA(interiorElementBoundaries),
                                                    IDATA(elementBoundaryElements),
                                                    IDATA(elementBoundaryLocalElementBoundaries),
                                                    DDATA(a),
                                                    DDATA(grad_phi),
                                                    DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryDiffusiveVelocity(PyObject* self, 
                                                            PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*a,*grad_phi,*velocity;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &a,
                       &grad_phi,
                       &velocity))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == ND(velocity));
      updateExteriorElementBoundaryDiffusiveVelocity(SHAPE(exteriorElementBoundaries)[0],
						     SHAPE(grad_phi)[1],
						     SHAPE(grad_phi)[2],
						     SHAPE(grad_phi)[3],
						     IDATA(exteriorElementBoundaries),
						     IDATA(elementBoundaryElements),
						     IDATA(elementBoundaryLocalElementBoundaries),
						     DDATA(a),
						     DDATA(grad_phi),
						     DDATA(velocity));
    }
  else
    {
      assert(ND(grad_phi) == ND(velocity));
      updateGlobalExteriorElementBoundaryDiffusiveVelocity(SHAPE(exteriorElementBoundaries)[0],
							   SHAPE(grad_phi)[1],
							   SHAPE(grad_phi)[2],
							   IDATA(exteriorElementBoundaries),
							   IDATA(elementBoundaryElements),
							   IDATA(elementBoundaryLocalElementBoundaries),
							   DDATA(a),
							   DDATA(grad_phi),
							   DDATA(velocity));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryDiffusiveVelocity_sd(PyObject* self, 
                                                            PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*a,*grad_phi,*velocity,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
                       &rowptr,
                       &colind,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &a,
                       &grad_phi,
                       &velocity))
    return NULL;
  if (ND(grad_phi) > 3)
    {
      assert(ND(grad_phi) == ND(velocity));
      updateExteriorElementBoundaryDiffusiveVelocity_sd(SHAPE(exteriorElementBoundaries)[0],
                                                        SHAPE(grad_phi)[1],
                                                        SHAPE(grad_phi)[2],
                                                        SHAPE(grad_phi)[3],
                                                        IDATA(rowptr),
                                                        IDATA(colind),
                                                        IDATA(exteriorElementBoundaries),
                                                        IDATA(elementBoundaryElements),
                                                        IDATA(elementBoundaryLocalElementBoundaries),
                                                        DDATA(a),
                                                        DDATA(grad_phi),
                                                        DDATA(velocity));
    }
  else
    {
      assert(ND(grad_phi) == ND(velocity));
      updateGlobalExteriorElementBoundaryDiffusiveVelocity_sd(SHAPE(exteriorElementBoundaries)[0],
                                                              SHAPE(grad_phi)[1],
                                                              SHAPE(grad_phi)[2],
                                                              IDATA(rowptr),
                                                              IDATA(colind),
                                                              IDATA(exteriorElementBoundaries),
                                                              IDATA(elementBoundaryElements),
                                                              IDATA(elementBoundaryLocalElementBoundaries),
                                                              DDATA(a),
                                                              DDATA(grad_phi),
                                                              DDATA(velocity));

    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateExteriorElementBoundaryAverageVelocity(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*velocity,*averageVelocity;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &velocity,
                       &averageVelocity))
    return NULL;
  calculateExteriorElementBoundaryAverageVelocity(SHAPE(exteriorElementBoundaries)[0],
                                                  SHAPE(velocity)[1],
                                                  SHAPE(velocity)[2],
                                                  SHAPE(velocity)[3],
                                                  IDATA(exteriorElementBoundaries),
                                                  IDATA(elementBoundaryElements),
                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                  DDATA(velocity),
                                                  DDATA(averageVelocity));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cfemIntegralsAccumulateExteriorElementPressureIntegrals(PyObject* self, 
							PyObject* args)
{
  PyObject *elementBoundaryMaterialTypes,*exteriorElementBoundaries,*pressure,*dS,*P,
    *boundaryMeasure;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &elementBoundaryMaterialTypes,
                       &exteriorElementBoundaries,
                       &pressure,
                       &dS,
                       &P,
		       &boundaryMeasure))
    return NULL;
  accumulateExteriorElementPressureIntegrals(SHAPE(exteriorElementBoundaries)[0],
					     SHAPE(pressure)[1],
					     IDATA(elementBoundaryMaterialTypes),
					     IDATA(exteriorElementBoundaries),
					     DDATA(pressure),
					     DDATA(dS),
					     DDATA(P),
					     DDATA(boundaryMeasure));
  Py_INCREF(Py_None);
  return Py_None;
}

/*mwf added temporarily*/
static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryShockCapturingVelocity(PyObject* self, 
								 PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*numDiff,
    *grad_u,*velocity;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &numDiff,
                       &grad_u,
                       &velocity))
    return NULL;
  updateInteriorElementBoundaryShockCapturingVelocity(SHAPE(interiorElementBoundaries)[0],
						      SHAPE(grad_u)[1],
						      SHAPE(grad_u)[2],
						      SHAPE(numDiff)[1],
						      SHAPE(grad_u)[3],
						      IDATA(interiorElementBoundaries),
						      IDATA(elementBoundaryElements),
						      IDATA(elementBoundaryLocalElementBoundaries),
						      DDATA(numDiff),
						      DDATA(grad_u),
						      DDATA(velocity));
  Py_INCREF(Py_None);
  return Py_None;
}
/*mwf add temporarily*/
static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryShockCapturingVelocity(PyObject* self, 
								 PyObject* args)
{
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*numDiff,*grad_u,
    *velocity;
  if(!PyArg_ParseTuple(args,"OOOOOO",
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &numDiff,
                       &grad_u,
                       &velocity))
    return NULL;
  if (ND(grad_u) > 3)
    {
      assert(SHAPE(numDiff)[0] == SHAPE(grad_u)[0]);
      updateExteriorElementBoundaryShockCapturingVelocity(SHAPE(exteriorElementBoundaries)[0],
							  SHAPE(grad_u)[1],
							  SHAPE(grad_u)[2],
							  SHAPE(numDiff)[1],
							  SHAPE(grad_u)[3],
							  IDATA(exteriorElementBoundaries),
							  IDATA(elementBoundaryElements),
							  IDATA(elementBoundaryLocalElementBoundaries),
							  DDATA(numDiff),
							  DDATA(grad_u),
							  DDATA(velocity));
    }
  else
    {
      assert(SHAPE(numDiff)[0] == SHAPE(grad_u)[0]);
      updateGlobalExteriorElementBoundaryShockCapturingVelocity(SHAPE(exteriorElementBoundaries)[0],
								SHAPE(grad_u)[1],
								SHAPE(grad_u)[2],
								IDATA(exteriorElementBoundaries),
								IDATA(elementBoundaryElements),
								IDATA(elementBoundaryLocalElementBoundaries),
								DDATA(numDiff),
								DDATA(grad_u),
								DDATA(velocity));
 
    }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateInteriorElementBoundaryAverageVelocity(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*velocity,*averageVelocity;
  if(!PyArg_ParseTuple(args,"OOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &velocity,
                       &averageVelocity))
    return NULL;
  calculateInteriorElementBoundaryAverageVelocity(SHAPE(interiorElementBoundaries)[0],
                                                  SHAPE(velocity)[1],
                                                  SHAPE(velocity)[2],
                                                  SHAPE(velocity)[3],
                                                  IDATA(interiorElementBoundaries),
                                                  IDATA(elementBoundaryElements),
                                                  IDATA(elementBoundaryLocalElementBoundaries),
                                                  DDATA(velocity),
                                                  DDATA(averageVelocity));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateConservationResidual(PyObject* self,
					   PyObject* args)
{
  
  PyObject *n,*dS_u,*elementResidual,*velocity,
    *conservationResidual;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &n,
		       &dS_u,
		       &elementResidual,
		       &velocity,
		       &conservationResidual))
    return NULL;

  calculateConservationResidual(SHAPE(elementResidual)[0],
				SHAPE(elementResidual)[1],
				SHAPE(n)[1],
				SHAPE(n)[2],
				SHAPE(n)[3],
				DDATA(n),
				DDATA(dS_u),
				DDATA(elementResidual),
				DDATA(velocity),
				DDATA(conservationResidual));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateConservationResidualDG(PyObject* self,
                                             PyObject* args)
{
  
  PyObject *elementResidual,*conservationResidual;
  if(!PyArg_ParseTuple(args,"OO",
		       &elementResidual,
		       &conservationResidual))
    return NULL;

  calculateConservationResidualDG(SHAPE(elementResidual)[0],
                                  SHAPE(elementResidual)[1],
                                  DDATA(elementResidual),
                                  DDATA(conservationResidual));
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCopyGlobalElementBoundaryVelocityToElementBoundary(PyObject* self, 
							   PyObject* args)
{
  PyObject* interiorElementBoundaries;
  PyObject* exteriorElementBoundaries;
  PyObject* elementBoundaryElements;
  PyObject* elementBoundaryLocalElementBoundaries;
  PyObject* ebq_global_v;
  PyObject* ebq_v;

  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &interiorElementBoundaries,
		       &exteriorElementBoundaries,
		       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebq_global_v,
		       &ebq_v))

    return NULL;
  copyGlobalElementBoundaryVelocityToElementBoundary(SHAPE(ebq_v)[0],
						     SHAPE(interiorElementBoundaries)[0],
						     SHAPE(exteriorElementBoundaries)[0],
						     SHAPE(ebq_global_v)[0],
						     SHAPE(ebq_v)[1],
						     SHAPE(ebq_v)[2],
						     SHAPE(ebq_v)[3],
						     IDATA(interiorElementBoundaries),
						     IDATA(exteriorElementBoundaries),
						     IDATA(elementBoundaryElements),
						     IDATA(elementBoundaryLocalElementBoundaries),
						     DDATA(ebq_global_v),
						     DDATA(ebq_v));


  Py_INCREF(Py_None); 
  return Py_None;
}
static PyObject*
cfemIntegralsLoadBoundaryFluxIntoGlobalElementBoundaryVelocity(PyObject* self, 
							   PyObject* args)
{
  PyObject* exteriorElementBoundaries;
  PyObject* fluxElementBoundaries;
  PyObject* n;
  PyObject* flux;
  PyObject* velocity;
  double updateCoef;
  if(!PyArg_ParseTuple(args,"OOOOdO",
		       &exteriorElementBoundaries,
		       &fluxElementBoundaries,
		       &n,
		       &flux,
		       &updateCoef,
		       &velocity))


    return NULL;
  loadBoundaryFluxIntoGlobalElementBoundaryVelocity(SHAPE(exteriorElementBoundaries)[0],
						    SHAPE(n)[1],
						    SHAPE(n)[2],
						    IDATA(exteriorElementBoundaries),
						    IDATA(fluxElementBoundaries),
						    DDATA(n),
						    DDATA(flux),
						    updateCoef,
						    DDATA(velocity));


  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateInteriorNumericalTrace_Potential(PyObject* self,
                                                       PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *phi,
    *dphi,
    *phi_trace,
    *dphi_trace_left,
    *dphi_trace_right;
  if(!PyArg_ParseTuple(args,"OOOOOOOO",
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &phi,
                       &dphi,
                       &phi_trace,
                       &dphi_trace_left,
                       &dphi_trace_right))
    return NULL;
  calculateInteriorNumericalTrace_Potential(SHAPE(interiorElementBoundaries)[0],
                                            SHAPE(phi)[1],
                                            SHAPE(phi)[2],
                                            IDATA(interiorElementBoundaries),
                                            IDATA(elementBoundaryElements),
                                            IDATA(elementBoundaryLocalElementBoundaries),
                                            DDATA(phi),
                                            DDATA(dphi),
                                            DDATA(phi_trace),
                                            DDATA(dphi_trace_left),
                                            DDATA(dphi_trace_right));
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateExteriorNumericalTrace_Potential(PyObject* self,
                                                       PyObject* args)
{ 
  PyObject *isDOFBoundary,
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *phi_bc,
    *phi,
    *dphi,
    *phi_trace,
    *dphi_trace_left;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &isDOFBoundary,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &phi_bc,
                       &phi,
                       &dphi,
                       &phi_trace,
                       &dphi_trace_left))
    return NULL;
  calculateExteriorNumericalTrace_Potential(IDATA(isDOFBoundary),
                                            SHAPE(exteriorElementBoundaries)[0],
                                            SHAPE(phi)[1],
                                            SHAPE(phi)[2],
                                            IDATA(exteriorElementBoundaries),
                                            IDATA(elementBoundaryElements),
                                            IDATA(elementBoundaryLocalElementBoundaries),
                                            DDATA(phi_bc),
                                            DDATA(phi),
                                            DDATA(dphi),
                                            DDATA(phi_trace),
                                            DDATA(dphi_trace_left));
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundary_MixedForm_weak(PyObject* self,
                                                          PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *phi_trace,
    *w_dS,
    *b;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &phi_trace,
                       &w_dS,
                       &b))
    return NULL;
  updateInteriorElementBoundary_MixedForm_weak(SHAPE(interiorElementBoundaries)[0],
                                               SHAPE(w_dS)[1],
                                               SHAPE(w_dS)[2],
                                               SHAPE(w_dS)[3],
                                               SHAPE(n)[3],
                                               IDATA(interiorElementBoundaries),
                                               IDATA(elementBoundaryElements),
                                               IDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(phi_trace),
                                               DDATA(w_dS),
                                               DDATA(b)); 
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundary_MixedForm_weakJacobian(PyObject* self,
                                                                   PyObject* args)
{ 
  PyObject *interiorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *dphi_trace_left,
    *dphi_trace_right,
    *v,
    *w_dS,
    *db,
    *db_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &dphi_trace_left,
                       &dphi_trace_right,
                       &v,
                       &w_dS,
                       &db,
                       &db_eb))
    return NULL;
  updateInteriorElementBoundary_MixedForm_weakJacobian(SHAPE(interiorElementBoundaries)[0],
                                                        SHAPE(v)[1],
                                                        SHAPE(v)[2],
                                                        SHAPE(v)[3],
                                                        SHAPE(n)[3],
                                                        IDATA(interiorElementBoundaries),
                                                        IDATA(elementBoundaryElements),
                                                        IDATA(elementBoundaryLocalElementBoundaries),
                                                        DDATA(n),
                                                        DDATA(dphi_trace_left),
                                                        DDATA(dphi_trace_right),
                                                        DDATA(v),
                                                        DDATA(w_dS),
                                                        DDATA(db),
                                                        DDATA(db_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundary_MixedForm_weak(PyObject* self,
                                                          PyObject* args)
{ 
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *phi_trace,
    *w_dS,
    *b;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &phi_trace,
                       &w_dS,
                       &b))
    return NULL;
  updateExteriorElementBoundary_MixedForm_weak(SHAPE(exteriorElementBoundaries)[0],
                                               SHAPE(w_dS)[1],
                                               SHAPE(w_dS)[2],
                                               SHAPE(w_dS)[3],
                                               SHAPE(n)[3],
                                               IDATA(exteriorElementBoundaries),
                                               IDATA(elementBoundaryElements),
                                               IDATA(elementBoundaryLocalElementBoundaries),
                                               DDATA(n),
                                               DDATA(phi_trace),
                                               DDATA(w_dS),
                                               DDATA(b)); 
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundary_MixedForm_weakJacobian(PyObject* self,
                                                                   PyObject* args)
{ 
  PyObject *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *n,
    *dphi_trace_left,
    *v,
    *w_dS,
    *db,
    *db_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &n,
                       &dphi_trace_left,
                       &v,
                       &w_dS,
                       &db,
                       &db_eb))
    return NULL;
  updateExteriorElementBoundary_MixedForm_weakJacobian(SHAPE(exteriorElementBoundaries)[0],
                                                       SHAPE(v)[1],
                                                       SHAPE(v)[2],
                                                       SHAPE(v)[3],
                                                       SHAPE(n)[3],
                                                       IDATA(exteriorElementBoundaries),
                                                       IDATA(elementBoundaryElements),
                                                       IDATA(elementBoundaryLocalElementBoundaries),
                                                       DDATA(n),
                                                       DDATA(dphi_trace_left),
                                                       DDATA(v),
                                                       DDATA(w_dS),
                                                       DDATA(db),
                                                       DDATA(db_eb));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cfemIntegralsUpdatePotential_MixedForm_weak(PyObject* self,
                                            PyObject* args)
{ 
  PyObject *phi,
    *grad_w_dV,
    *b;
  if(!PyArg_ParseTuple(args,"OOO",
		       &phi,
                       &grad_w_dV,
                       &b))
    return NULL;
  updatePotential_MixedForm_weak(SHAPE(grad_w_dV)[0],
                                 SHAPE(grad_w_dV)[1],
                                 SHAPE(grad_w_dV)[2],
                                 SHAPE(grad_w_dV)[3],
                                 DDATA(phi),
                                 DDATA(grad_w_dV),
                                 DDATA(b));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdatePotential_MixedForm_weak_gwvd(PyObject* self,
                                            PyObject* args)
{ 
  double epsilon;
  PyObject *phi, *w_dV, *mf,*grad_w_dV,*b;
  if(!PyArg_ParseTuple(args,"OOOdOO",
		       &phi,
                       &grad_w_dV,
		       &w_dV,
		       &epsilon,
		       &mf,
                       &b))
    return NULL;
  updatePotential_MixedForm_weak_gwvd(SHAPE(grad_w_dV)[0],
                                 SHAPE(grad_w_dV)[1],
                                 SHAPE(grad_w_dV)[2],
                                 SHAPE(grad_w_dV)[3],
				 epsilon,
                                 DDATA(phi),
				 DDATA(w_dV),
                                 DDATA(grad_w_dV),
				 DDATA(b),   
				 DDATA(mf));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdatePotential_MixedForm_weakJacobian(PyObject* self,
                                                    PyObject* args)
{ 
  PyObject *dphi,
    *v,
    *grad_w_dV,
    *b;
  if(!PyArg_ParseTuple(args,"OOOO",
		       &dphi,
		       &v,
                       &grad_w_dV,
                       &b))
    return NULL;
  updatePotential_MixedForm_weakJacobian(SHAPE(grad_w_dV)[0],
                                         SHAPE(grad_w_dV)[1],
                                         SHAPE(grad_w_dV)[2],
                                         SHAPE(grad_w_dV)[3],
                                         DDATA(dphi),
                                         DDATA(v),
                                         DDATA(grad_w_dV),
                                         DDATA(b));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateVelocityQuadrature_MixedForm(PyObject* self,
                                                   PyObject* args)
{ 
  PyObject *A_inv,
    *b,
    *v,
    *V,
    *qv,
    *qV;
  if(!PyArg_ParseTuple(args,"OOOOOO",
		       &A_inv,
                       &b,
                       &v,
                       &V,
                       &qv,
                       &qV))
    return NULL;
  calculateVelocityQuadrature_MixedForm(SHAPE(v)[0],
                                        SHAPE(v)[1],
                                        SHAPE(v)[2],
                                        SHAPE(v)[3],
                                        SHAPE(b)[1],
                                        SHAPE(qv)[1],
                                        DDATA(A_inv),
                                        DDATA(b),
                                        DDATA(v),
                                        DDATA(V),
                                        DDATA(qv),
                                        DDATA(qV));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cfemIntegralsCalculateVelocityQuadrature_MixedForm2(PyObject* self,
                                                   PyObject* args)
{ 
  PyObject *qa,
    *qw_dV,
    *b,
    *v,
    *V,
    *qv,
    *qV;
  if(!PyArg_ParseTuple(args,"OOOOOOO",
		       &qa,
		       &qw_dV,
                       &b,
                       &v,
                       &V,
                       &qv,
                       &qV))
    return NULL;
      calculateVelocityQuadrature_MixedForm2(SHAPE(v)[0],
					     SHAPE(v)[1],
					     SHAPE(v)[2],
					     SHAPE(v)[3],
					     SHAPE(b)[1],
					     SHAPE(qv)[1],
					     DDATA(qa),
					     DDATA(qw_dV),
					     DDATA(b),
					     DDATA(v),
					     DDATA(V),
					     DDATA(qv),
					     DDATA(qV));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject*
cfemIntegralsCalculateVelocityQuadrature_MixedForm2_sd(PyObject* self,
                                                   PyObject* args)
{ 
  PyObject *qa,
    *qw_dV,
    *b,
    *v,
    *V,
    *qv,
    *qV,
    *rowptr,
    *colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &rowptr,
		       &colind,
		       &qa,
		       &qw_dV,
                       &b,
                       &v,
                       &V,
                       &qv,
                       &qV))
    return NULL;
      calculateVelocityQuadrature_MixedForm2_sd(SHAPE(v)[0],
						SHAPE(v)[1],
						SHAPE(v)[2],
						SHAPE(v)[3],
						SHAPE(b)[1],
						SHAPE(qv)[1],
						IDATA(rowptr),
						IDATA(colind),
						DDATA(qa),
						DDATA(qw_dV),
						DDATA(b),
						DDATA(v),
						DDATA(V),
						DDATA(qv),
						DDATA(qV));
      Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cgwvdNumericalFluxCalculateVelocityQuadrature_MixedForm2_vdof_sd(PyObject* self,
                                                   PyObject* args)
{ 
  PyObject *qa,
    *qw_dV,
    *b,
    *v,
    *V,
    *qv,
    *qV,
    *rowptr,
    *colind,
    *vel_dofs;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &qa,
		       &qw_dV,
                       &b,
                       &v,
                       &V,
                       &qv,
                       &qV,
		       &vel_dofs))
    return NULL;
      calculateVelocityQuadrature_MixedForm2_vdof_sd(SHAPE(v)[0],
						SHAPE(v)[1],
						SHAPE(v)[2],
						SHAPE(v)[3],
						SHAPE(b)[1],
						SHAPE(qv)[1],
						IDATA(rowptr),
						IDATA(colind),
						DDATA(qa),
						DDATA(qw_dV),
						DDATA(b),
						DDATA(v),
						DDATA(V),
						DDATA(qv),
						DDATA(qV),
						DDATA(vel_dofs));
      Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateVelocityQuadrature_MixedForm_Jacobian(PyObject* self,
                                                            PyObject* args)
{ 
  PyObject *A_inv,
    *db,
    *db_eb,
    *v,
    *DV,
    *DV_eb,
    *qv,
    *qDV,
    *qDV_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
		       &A_inv,
                       &db,
                       &db_eb,
                       &v,
                       &DV,
                       &DV_eb,
                       &qv,
                       &qDV,
                       &qDV_eb))
    return NULL;
  calculateVelocityQuadrature_MixedForm_Jacobian(SHAPE(v)[0],
						  SHAPE(v)[1],
						  SHAPE(v)[2],
						  SHAPE(v)[3],
						  SHAPE(db)[1],
						  SHAPE(qv)[1],
						  DDATA(A_inv),
						  DDATA(db),
						  DDATA(db_eb),
						  DDATA(v),
						  DDATA(DV),
						  DDATA(DV_eb),
						  DDATA(qv),
						  DDATA(qDV),
						  DDATA(qDV_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateVelocityQuadrature_MixedForm2_Jacobian(PyObject* self,
                                                            PyObject* args)
{ 
  PyObject *qa,
    *qw_dV,
    *db,
    *db_eb,
    *v,
    *DV,
    *DV_eb,
    *qv,
    *qDV,
    *qDV_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
		       &qa,
		       &qw_dV,
                       &db,
                       &db_eb,
                       &v,
                       &DV,
                       &DV_eb,
                       &qv,
                       &qDV,
                       &qDV_eb))
    return NULL;
  calculateVelocityQuadrature_MixedForm2_Jacobian(SHAPE(v)[0],
						  SHAPE(v)[1],
						  SHAPE(v)[2],
						  SHAPE(v)[3],
						  SHAPE(db)[1],
						  SHAPE(qv)[1],
						  DDATA(qa),
						  DDATA(qw_dV),
						  DDATA(db),
						  DDATA(db_eb),
						  DDATA(v),
						  DDATA(DV),
						  DDATA(DV_eb),
						  DDATA(qv),
						  DDATA(qDV),
						  DDATA(qDV_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateVelocityQuadrature_MixedForm2_Jacobian_sd(PyObject* self,
								PyObject* args)
{ 
  PyObject *qa,
    *qw_dV,
    *db,
    *db_eb,
    *v,
    *DV,
    *DV_eb,
    *qv,
    *qDV,
    *qDV_eb,
    *rowptr,
    *colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOO",
		       &rowptr,
		       &colind,
		       &qa,
		       &qw_dV,
                       &db,
                       &db_eb,
                       &v,
                       &DV,
                       &DV_eb,
                       &qv,
                       &qDV,
                       &qDV_eb))
    return NULL;
  calculateVelocityQuadrature_MixedForm2_Jacobian_sd(SHAPE(v)[0],
						     SHAPE(v)[1],
						     SHAPE(v)[2],
						     SHAPE(v)[3],
						     SHAPE(db)[1],
						     SHAPE(qv)[1],
						     IDATA(rowptr),
						     IDATA(colind),
						     DDATA(qa),
						     DDATA(qw_dV),
						     DDATA(db),
						     DDATA(db_eb),
						     DDATA(v),
						     DDATA(DV),
						     DDATA(DV_eb),
						     DDATA(qv),
						     DDATA(qDV),
						     DDATA(qDV_eb));
  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
cfemIntegralsCalculateVelocityProjectionMatrixLDG(PyObject* self,
                                                  PyObject* args)
{
  PyObject *vXw_dV,
    *A_inv;
  if(!PyArg_ParseTuple(args,"OO",
		       &vXw_dV,
                       &A_inv))
    return NULL;
  calculateVelocityProjectionMatrixLDG(SHAPE(vXw_dV)[0],
                                       SHAPE(vXw_dV)[1],
                                       SHAPE(vXw_dV)[2],
                                       DDATA(vXw_dV),
                                       DDATA(A_inv));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_MixedForm_weak(PyObject* self, 
                                                             PyObject* args)
{
  PyObject *a,*qV,*grad_w_dV,*residual;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &a,
                       &qV,
                       &grad_w_dV,
                       &residual))
    return NULL;
  updateDiffusion_MixedForm_weak(SHAPE(grad_w_dV)[0],
                       SHAPE(grad_w_dV)[1],
                       SHAPE(grad_w_dV)[2],
                       SHAPE(grad_w_dV)[3],
                       DDATA(a),
                       DDATA(qV),
                       DDATA(grad_w_dV),
                       DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusion_MixedForm_weak_sd(PyObject* self, 
                                                                PyObject* args)
{
  PyObject *a,*qV,*grad_w_dV,*residual,*rowptr,*colind,*velocity;
  int rho_split=0;
  
  if(!PyArg_ParseTuple(args,"OOOOOOO|i",
		       &rowptr,
                       &colind,
                       &a,
                       &qV,
                       &grad_w_dV,
		       &velocity,
                       &residual,
		       &rho_split))
    return NULL;
  updateDiffusion_MixedForm_weak_sd(SHAPE(grad_w_dV)[0],
                                    SHAPE(grad_w_dV)[1],
                                    SHAPE(grad_w_dV)[2],
                                    SHAPE(grad_w_dV)[3],
				    rho_split,
                                    IDATA(rowptr),
                                    IDATA(colind),
                                    DDATA(a),
                                    DDATA(qV),
                                    DDATA(grad_w_dV),
				    DDATA(velocity),
                                    DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_MixedForm_weak(PyObject* self, 
                                                                     PyObject* args)
{
  PyObject *a,*da,*v,*qV,*qDV,*qDV_eb,*grad_w_dV,*jacobian,*jacobian_eb;
  if(!PyArg_ParseTuple(args,"OOOOOOOOO",
                       &a,
                       &da,
                       &qV,
                       &qDV,
                       &qDV_eb,
                       &grad_w_dV,
                       &v,
                       &jacobian,
                       &jacobian_eb))
    return NULL;
  updateDiffusionJacobian_MixedForm_weak(SHAPE(qDV_eb)[0],
                                         SHAPE(qDV_eb)[1],
                                         SHAPE(qDV_eb)[2],
                                         SHAPE(qDV_eb)[3],
                                         SHAPE(grad_w_dV)[2],
                                         SHAPE(grad_w_dV)[3],
                                         DDATA(a),
                                         DDATA(da),
                                         DDATA(qV),
                                         DDATA(qDV),
                                         DDATA(qDV_eb),
                                         DDATA(grad_w_dV),
                                         DDATA(v),
                                         DDATA(jacobian),
                                         DDATA(jacobian_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateDiffusionJacobian_MixedForm_weak_sd(PyObject* self, 
                                                                     PyObject* args)
{
  PyObject *a,*da,*v,*qV,*qDV,*qDV_eb,*grad_w_dV,*jacobian,*jacobian_eb,*rowptr,*colind;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
                       &rowptr,
                       &colind,
                       &a,
                       &da,
                       &qV,
                       &qDV,
                       &qDV_eb,
                       &grad_w_dV,
                       &v,
                       &jacobian,
                       &jacobian_eb))
    return NULL;
  updateDiffusionJacobian_MixedForm_weak_sd(SHAPE(qDV_eb)[0],
                                            SHAPE(qDV_eb)[1],
                                            SHAPE(qDV_eb)[2],
                                            SHAPE(qDV_eb)[3],
                                            SHAPE(grad_w_dV)[2],
                                            SHAPE(grad_w_dV)[3],
                                            IDATA(rowptr),
                                            IDATA(colind),
                                            DDATA(a),
                                            DDATA(da),
                                            DDATA(qV),
                                            DDATA(qDV),
                                            DDATA(qDV_eb),
                                            DDATA(grad_w_dV),
                                            DDATA(v),
                                            DDATA(jacobian),
                                            DDATA(jacobian_eb));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsEstimate_mt(PyObject* self,
                         PyObject* args)
{
  PyObject *v,
    *vXw_dV,
    *elementSpatialResidual,
    *mt;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &v,
		       &vXw_dV,
                       &elementSpatialResidual,
                       &mt))
    return NULL;
  estimate_mt(SHAPE(vXw_dV)[0],
              SHAPE(vXw_dV)[1],
              SHAPE(vXw_dV)[2],
              DDATA(v),
              DDATA(vXw_dV),
              DDATA(elementSpatialResidual),
              DDATA(mt));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsEstimate_mt_lowmem(PyObject* self,
                                PyObject* args)
{
  PyObject *v,
    *w_dV,
    *elementSpatialResidual,
    *mt;
  if(!PyArg_ParseTuple(args,"OOOO",
                       &v,
		       &w_dV,
                       &elementSpatialResidual,
                       &mt))
    return NULL;
  estimate_mt_lowmem(SHAPE(w_dV)[0],
                     SHAPE(w_dV)[1],
                     SHAPE(w_dV)[2],
                     DDATA(v),
                     DDATA(w_dV),
                     DDATA(elementSpatialResidual),
                     DDATA(mt));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCopyExteriorElementBoundaryValuesFromElementBoundaryValues(PyObject* self,
									PyObject* args)
{
  PyObject *exteriorElementBoundaries, *elementBoundaryElements, *elementBoundaryLocalElementBoundaries,
    *ebq_val,*ebqe_val;
  int nValuesPerQuadraturePoint = 1,i;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebq_val,
		       &ebqe_val))
    return NULL;
  if (ND(ebq_val) > 3)
    for (i = 3; i < ND(ebq_val); i++)
      nValuesPerQuadraturePoint *= SHAPE(ebq_val)[i];


  copyExteriorElementBoundaryValuesFromElementBoundaryValues(SHAPE(exteriorElementBoundaries)[0],
							     SHAPE(ebq_val)[0],
							     SHAPE(ebq_val)[1],
							     SHAPE(ebq_val)[2],
							     nValuesPerQuadraturePoint,
							     IDATA(exteriorElementBoundaries),
							     IDATA(elementBoundaryElements),
							     IDATA(elementBoundaryLocalElementBoundaries),
							     DDATA(ebq_val),
							     DDATA(ebqe_val));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCopyExteriorElementBoundaryValuesToElementBoundaryValues(PyObject* self,
									PyObject* args)
{
  PyObject *exteriorElementBoundaries, *elementBoundaryElements, *elementBoundaryLocalElementBoundaries,
    *ebq_val,*ebqe_val;
  int nValuesPerQuadraturePoint = 1,i;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebqe_val,
		       &ebq_val))
    return NULL;
  if (ND(ebq_val) > 3)
    for (i = 3; i < ND(ebq_val); i++)
      nValuesPerQuadraturePoint *= SHAPE(ebq_val)[i];


  copyExteriorElementBoundaryValuesToElementBoundaryValues(SHAPE(exteriorElementBoundaries)[0],
							   SHAPE(ebq_val)[0],
							   SHAPE(ebq_val)[1],
							   SHAPE(ebq_val)[2],
							   nValuesPerQuadraturePoint,
							   IDATA(exteriorElementBoundaries),
							   IDATA(elementBoundaryElements),
							   IDATA(elementBoundaryLocalElementBoundaries),
							   DDATA(ebqe_val),
							   DDATA(ebq_val));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCopyExteriorElementBoundaryValuesToGlobalElementBoundaryValues(PyObject* self,
									    PyObject* args)
{
  PyObject *exteriorElementBoundaries, *elementBoundaryElements, *elementBoundaryLocalElementBoundaries,
    *ebq_global_val,*ebqe_val;
  int nValuesPerQuadraturePoint = 1,i;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebqe_val,
		       &ebq_global_val))
    return NULL;
  if (ND(ebq_global_val) > 2)
    for (i = 2; i < ND(ebq_global_val); i++)
      nValuesPerQuadraturePoint *= SHAPE(ebq_global_val)[i];


  copyExteriorElementBoundaryValuesToGlobalElementBoundaryValues(SHAPE(exteriorElementBoundaries)[0],
								 SHAPE(ebq_global_val)[1],
								 nValuesPerQuadraturePoint,
								 IDATA(exteriorElementBoundaries),
								 IDATA(elementBoundaryElements),
								 IDATA(elementBoundaryLocalElementBoundaries),
								 DDATA(ebqe_val),
								 DDATA(ebq_global_val));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsCopyExteriorElementBoundaryValuesFromGlobalElementBoundaryValues(PyObject* self,
									      PyObject* args)
{
  PyObject *exteriorElementBoundaries, *elementBoundaryElements, *elementBoundaryLocalElementBoundaries,
    *ebq_global_val,*ebqe_val;
  int nValuesPerQuadraturePoint = 1,i;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebq_global_val,
		       &ebqe_val))
    return NULL;
  if (ND(ebq_global_val) > 2)
    for (i = 2; i < ND(ebq_global_val); i++)
      nValuesPerQuadraturePoint *= SHAPE(ebq_global_val)[i];

  copyExteriorElementBoundaryValuesFromGlobalElementBoundaryValues(SHAPE(exteriorElementBoundaries)[0],
								   SHAPE(ebq_global_val)[1],
								   nValuesPerQuadraturePoint,
								   IDATA(exteriorElementBoundaries),
								   IDATA(elementBoundaryElements),
								   IDATA(elementBoundaryLocalElementBoundaries),
								   DDATA(ebq_global_val),
								   DDATA(ebqe_val));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsScalarDomainIntegral(PyObject* self,
                                  PyObject* args)
{
  double integral=0.0;
  int nElements;
  PyObject *dV,*nValueArray;
  if(!PyArg_ParseTuple(args,"OOi",
                       &dV,
                       &nValueArray,
		       &nElements))
    return NULL;
  integral = scalarDomainIntegral(nElements,
				  SHAPE(dV)[1],
				  DDATA(dV),
				  DDATA(nValueArray));
  return Py_BuildValue("d",integral);
}

static PyObject*
cfemIntegralsScalarHeavisideDomainIntegral(PyObject* self,
					   PyObject* args)
{
  int nElements;
  double integral=0.0;
  PyObject *dV,*nValueArray;
  if(!PyArg_ParseTuple(args,"OOi",
                       &dV,
                       &nValueArray,
		       &nElements))
    return NULL;
  integral = scalarHeavisideDomainIntegral(nElements,
					   SHAPE(dV)[1],
					   DDATA(dV),
					   DDATA(nValueArray));
  return Py_BuildValue("d",integral);
}

static PyObject*
cfemIntegralsScalarSmoothedHeavisideDomainIntegral(PyObject* self,
						   PyObject* args)
{
  int nElements;
  double integral=0.0,epsFact;
  PyObject *dV,*nValueArray,*elementDiameter;
  if(!PyArg_ParseTuple(args,"dOOOi",
		       &epsFact,
		       &elementDiameter,
                       &dV,
                       &nValueArray,
		       &nElements))
    return NULL;
  integral = scalarSmoothedHeavisideDomainIntegral(nElements,
						   SHAPE(dV)[1],
						   epsFact,
						   DDATA(elementDiameter),
						   DDATA(dV),
						   DDATA(nValueArray));
  return Py_BuildValue("d",integral);
}

static PyObject*
cfemIntegralsFluxDomainBoundaryIntegral(PyObject* self,
                                        PyObject* args)
{
  int nElementBoundaries_owned;
  double integral=0.0;
  PyObject *flag,*dS,*nValueArray,*exteriorElementBoundariesArray;
  if(!PyArg_ParseTuple(args,"iOOOO",
		       &nElementBoundaries_owned,
                       &flag,
                       &exteriorElementBoundariesArray,
                       &dS,
                       &nValueArray))
    return NULL;
  integral = fluxDomainBoundaryIntegral(SHAPE(nValueArray)[0],
					nElementBoundaries_owned,
                                        SHAPE(nValueArray)[1],
                                        IDATA(flag),
                                        IDATA(exteriorElementBoundariesArray),
                                        DDATA(dS),
                                        DDATA(nValueArray));
  return Py_BuildValue("d",integral);
}

static PyObject*
cfemIntegralsFluxDomainBoundaryIntegralFromVector(PyObject* self,
						  PyObject* args)
{
  int nElementBoundaries_owned;
  double integral=0.0;
  PyObject *flag,*exteriorElementBoundaries,*dS,*nValueArray,*normal;
  if(!PyArg_ParseTuple(args,"iOOOOO",
		       &nElementBoundaries_owned,
		       &flag,
		       &exteriorElementBoundaries,
                       &dS,
                       &nValueArray,
		       &normal))
    return NULL;
  integral = fluxDomainBoundaryIntegralFromVector(SHAPE(normal)[0],
						  nElementBoundaries_owned,
						  SHAPE(normal)[1],
						  SHAPE(normal)[2],
						  IDATA(flag),
						  IDATA(exteriorElementBoundaries),
						  DDATA(dS),
						  DDATA(nValueArray),
						  DDATA(normal));
  return Py_BuildValue("d",integral);
}

static PyObject*
cfemIntegralsComputeC0P1InterpolantDGP12(PyObject* self,
					PyObject* args)
{
  PyObject *elementNodesArray,*nodeElementOffsets,*nodeElementsArray, *l2g, *dof, *nodalAverage;
  int dim_dof = 1;
  if(!PyArg_ParseTuple(args,"OOOOOO|i",
                       &elementNodesArray,
		       &nodeElementOffsets,
		       &nodeElementsArray,
		       &l2g,
                       &dof,
		       &nodalAverage,
		       &dim_dof))
    return NULL;

  computeC0P1InterpolantDGP12(SHAPE(elementNodesArray)[0],
			      SHAPE(nodeElementOffsets)[0]-1,
			      SHAPE(elementNodesArray)[1],
			      SHAPE(l2g)[1],
			      dim_dof,
			      IDATA(elementNodesArray),
			      IDATA(nodeElementOffsets),
			      IDATA(nodeElementsArray),
			      IDATA(l2g),
			      DDATA(dof),
			      DDATA(nodalAverage));
  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
cfemIntegralsComputeC0P1InterpolantDGP0(PyObject* self,
					PyObject* args)
{
  PyObject *elementNodesArray,*nodeElementOffsets,*nodeElementsArray, *l2g, *dof, *nodalAverage;
  int dim_dof = 1;
  if(!PyArg_ParseTuple(args,"OOOOOO|i",
                       &elementNodesArray,
		       &nodeElementOffsets,
		       &nodeElementsArray,
		       &l2g,
                       &dof,
		       &nodalAverage,
		       &dim_dof))
    return NULL;

  computeC0P1InterpolantDGP0(SHAPE(elementNodesArray)[0],
			     SHAPE(nodeElementOffsets)[0]-1,
			     SHAPE(elementNodesArray)[1],
			     SHAPE(l2g)[1],
			     dim_dof,
			     IDATA(elementNodesArray),
			     IDATA(nodeElementOffsets),
			     IDATA(nodeElementsArray),
			     IDATA(l2g),
			     DDATA(dof),
			     DDATA(nodalAverage));
  Py_INCREF(Py_None);
  return Py_None;

}
static PyObject*
cfemIntegralsComputeC0P1InterpolantNCP1(PyObject* self,
					PyObject* args)
{
  PyObject *elementNodesArray,*nodeElementOffsets,*nodeElementsArray, *l2g, *dof, *nodalAverage;
  int dim_dof = 1;
  if(!PyArg_ParseTuple(args,"OOOOOO|i",
                       &elementNodesArray,
		       &nodeElementOffsets,
		       &nodeElementsArray,
		       &l2g,
                       &dof,
		       &nodalAverage,
		       &dim_dof))
    return NULL;

  computeC0P1InterpolantNCP1(SHAPE(elementNodesArray)[0],
			     SHAPE(nodeElementOffsets)[0]-1,
			     SHAPE(elementNodesArray)[1],
			     SHAPE(l2g)[1],
			     dim_dof,
			     IDATA(elementNodesArray),
			     IDATA(nodeElementOffsets),
			     IDATA(nodeElementsArray),
			     IDATA(l2g),
			     DDATA(dof),
			     DDATA(nodalAverage));
  Py_INCREF(Py_None);
  return Py_None;

}


static PyObject*
cfemIntegralsCheckElementBoundaryAndExteriorElementBoundaryArraysSame(PyObject* self,
								      PyObject* args)
{
  double tolerance;
  PyObject *exteriorElementBoundaries, *elementBoundaryElements, *elementBoundaryLocalElementBoundaries,
    *ebq_val,*ebqe_val;
  int nValuesPerQuadraturePoint = 1,failed = 0,firstBadIndex= -1,i;
  if(!PyArg_ParseTuple(args,"dOOOOO",
		       &tolerance,
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebq_val,
		       &ebqe_val))
    return NULL;
  if (ND(ebqe_val) > 2)
    for (i=2; i < ND(ebqe_val); i++)
      nValuesPerQuadraturePoint *= SHAPE(ebqe_val)[i];


  failed = checkElementBoundaryAndExteriorElementBoundaryArraysSame(SHAPE(ebq_val)[1],
								    SHAPE(exteriorElementBoundaries)[0],
								    SHAPE(ebqe_val)[1],
								    nValuesPerQuadraturePoint,
								    tolerance,
								    IDATA(exteriorElementBoundaries),
								    IDATA(elementBoundaryElements),
								    IDATA(elementBoundaryLocalElementBoundaries),
								    DDATA(ebq_val),
								    DDATA(ebqe_val),
								    &firstBadIndex);

  return Py_BuildValue("(ii)",failed,firstBadIndex);
}
static PyObject*
cfemIntegralsCheckGlobalElementBoundaryAndExteriorElementBoundaryArraysSame(PyObject* self,
									    PyObject* args)
{
  double tolerance;
  PyObject *exteriorElementBoundaries, *elementBoundaryElements, *elementBoundaryLocalElementBoundaries,
    *ebq_global_val,*ebqe_val;
  int nValuesPerQuadraturePoint = 1,failed = 0,firstBadIndex= -1,i;
  if(!PyArg_ParseTuple(args,"dOOOOO",
		       &tolerance,
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
		       &elementBoundaryLocalElementBoundaries,
		       &ebq_global_val,
		       &ebqe_val))
    return NULL;
  if (ND(ebqe_val) > 2)
    for (i=2; i < ND(ebqe_val); i++)
      nValuesPerQuadraturePoint *= SHAPE(ebqe_val)[i];


  failed = checkGlobalElementBoundaryAndExteriorElementBoundaryArraysSame(SHAPE(exteriorElementBoundaries)[0],
									  SHAPE(ebqe_val)[1],
									  nValuesPerQuadraturePoint,
									  tolerance,
									  IDATA(exteriorElementBoundaries),
									  IDATA(elementBoundaryElements),
									  IDATA(elementBoundaryLocalElementBoundaries),
									  DDATA(ebq_global_val),
									  DDATA(ebqe_val),
									  &firstBadIndex);

  return Py_BuildValue("(ii)",failed,firstBadIndex);
}

static PyObject*
cfemIntegralsCalculateExteriorElementBoundaryStress3D(PyObject* self,
                                                      PyObject* args)
{
  PyObject *elementBoundaryMaterialTypes,
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *p,
    *mom_flux_vec_u,
    *mom_flux_vec_v,
    *mom_flux_vec_w,
    *dS,
    *n,
    *F;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOO",
                       &elementBoundaryMaterialTypes,
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &p,
                       &mom_flux_vec_u,
                       &mom_flux_vec_v,
                       &mom_flux_vec_w,
                       &dS,
                       &n,
                       &F))
    return NULL;
  calculateExteriorElementBoundaryStress3D(SHAPE(p)[0],
                                           SHAPE(p)[1],
                                           IDATA(elementBoundaryMaterialTypes),
                                           IDATA(exteriorElementBoundaries),
                                           IDATA(elementBoundaryElements),
                                           IDATA(elementBoundaryLocalElementBoundaries),
                                           DDATA(p),
                                           DDATA(mom_flux_vec_u),
                                           DDATA(mom_flux_vec_v),
                                           DDATA(mom_flux_vec_w),
                                           DDATA(dS),
                                           DDATA(n),
                                           DDATA(F));
  return Py_None;
}

static PyObject*
cfemIntegralsCalculateExteriorElementBoundaryStress2D(PyObject* self,
                                                      PyObject* args)
{
  PyObject *elementBoundaryMaterialTypes,
    *exteriorElementBoundaries,
    *elementBoundaryElements,
    *elementBoundaryLocalElementBoundaries,
    *p,
    *mom_flux_vec_u,
    *mom_flux_vec_v,
    *dS,
    *n,
    *F;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOO",
                       &elementBoundaryMaterialTypes,
		       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &p,
                       &mom_flux_vec_u,
                       &mom_flux_vec_v,
                       &dS,
                       &n,
                       &F))
    return NULL;
  calculateExteriorElementBoundaryStress2D(SHAPE(p)[0],
                                           SHAPE(p)[1],
                                           IDATA(elementBoundaryMaterialTypes),
                                           IDATA(exteriorElementBoundaries),
                                           IDATA(elementBoundaryElements),
                                           IDATA(elementBoundaryLocalElementBoundaries),
                                           DDATA(p),
                                           DDATA(mom_flux_vec_u),
                                           DDATA(mom_flux_vec_v),
                                           DDATA(dS),
                                           DDATA(n),
                                           DDATA(F));
  return Py_None;
}
static PyObject* cfemIntegralsCopyBetweenFreeUnknownsAndGlobalUnknowns(PyObject* self,PyObject* args)
{
  PyObject *globalDOFids,*freeDOFids,*free_u,*u;
  int copyFromFreeToGlobal,offset,stride;
  if (!PyArg_ParseTuple(args,
                        "iiiOOOO",
			&copyFromFreeToGlobal,
                        &offset,
			&stride,
                        &globalDOFids,
			&freeDOFids,
			&free_u,
			&u))
    return NULL;
  if (copyFromFreeToGlobal > 0)
    copyFreeUnknownsToGlobalUnknowns(SHAPE(globalDOFids)[0],
				     offset,
				     stride,
				     IDATA(globalDOFids),
				     IDATA(freeDOFids),
				     DDATA(free_u),
				     DDATA(u));
  else if  (copyFromFreeToGlobal == 0)
    copyGlobalUnknownsToFreeUnknowns(SHAPE(globalDOFids)[0],
				     offset,
				     stride,
				     IDATA(globalDOFids),
				     IDATA(freeDOFids),
				     DDATA(u),
				     DDATA(free_u));
  else
    {
      printf("error copyFromFreeToGlobal = %d not recognized quitting\n",copyFromFreeToGlobal);
      exit(1);
    }
  Py_INCREF(Py_None); 
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryDiffusionAdjoint(PyObject* self, 
							   PyObject* args)
{
  double sigma;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*u,*n,*a,*grad_w,*dS,*residual;
  if(!PyArg_ParseTuple(args,"OOOdOOOOOO",
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &sigma,
		       &u,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
		       &residual))
    return NULL;
  updateInteriorElementBoundaryDiffusionAdjoint(SHAPE(interiorElementBoundaries)[0],
						SHAPE(grad_w)[1],
						SHAPE(grad_w)[2],
						SHAPE(grad_w)[3],
						SHAPE(grad_w)[4],
						IDATA(interiorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						sigma,
						DDATA(u),
						DDATA(n),
						DDATA(a),
						DDATA(grad_w),
						DDATA(dS),
						DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryDiffusionAdjoint(PyObject* self, 
							   PyObject* args)
{
  double sigma;
  PyObject *isDOFBoundary,*exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*u,*ub,*n,*a,*grad_w,*dS,*residual;
  if(!PyArg_ParseTuple(args,"OOOOdOOOOOOO",
		       &isDOFBoundary,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &sigma,
		       &u,
		       &ub,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
		       &residual))
    return NULL;
  updateExteriorElementBoundaryDiffusionAdjoint(SHAPE(exteriorElementBoundaries)[0],
						SHAPE(grad_w)[1],
						SHAPE(grad_w)[2],
						SHAPE(grad_w)[3],
						IDATA(isDOFBoundary),
						IDATA(exteriorElementBoundaries),
						IDATA(elementBoundaryElements),
						IDATA(elementBoundaryLocalElementBoundaries),
						sigma,
						DDATA(u),
						DDATA(ub),
						DDATA(n),
						DDATA(a),
						DDATA(grad_w),
						DDATA(dS),
						DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense(PyObject* self, 
										   PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  double sigma;
  PyObject *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *v,*n,*a,*grad_w,*dS,*jac;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOdOOOOOO",
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &sigma,
		       &v,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense(SHAPE(interiorElementBoundaries)[0],
									SHAPE(grad_w)[1],
									SHAPE(grad_w)[2],
									SHAPE(grad_w)[3],
									SHAPE(v)[3],
									SHAPE(n)[3],
									offset_r,
									stride_r,
									offset_u,
									stride_u,
									nFreeVDOF_global,
									IDATA(interiorElementBoundaries),
									IDATA(elementBoundaryElements),
									IDATA(elementBoundaryLocalElementBoundaries),
									IDATA(nFreeDOF_element_r),
									IDATA(nFreeDOF_element_u),
									IDATA(freeLocal_r),
									IDATA(freeGlobal_r),
									IDATA(freeLocal_u),
									IDATA(freeGlobal_u),
									sigma,
									DDATA(v),
									DDATA(n),
									DDATA(a),
									DDATA(grad_w),
									DDATA(dS),
									DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense(PyObject* self, 
										   PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  double sigma;
  PyObject *exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *isDOFBoundary,*v,*n,*a,*grad_w,*dS,*jac;
  if(!PyArg_ParseTuple(args,"iiiiiOOOOOOOOOOdOOOOOO",
                       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &isDOFBoundary,
		       &sigma,
		       &v,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense(SHAPE(exteriorElementBoundaries)[0],
									SHAPE(grad_w)[1],
									SHAPE(grad_w)[2],
									SHAPE(v)[2],
									SHAPE(n)[2],
									offset_r,
									stride_r,
									offset_u,
									stride_u,
									nFreeVDOF_global,
									IDATA(exteriorElementBoundaries),
									IDATA(elementBoundaryElements),
									IDATA(elementBoundaryLocalElementBoundaries),
									IDATA(nFreeDOF_element_r),
									IDATA(nFreeDOF_element_u),
									IDATA(freeLocal_r),
									IDATA(freeGlobal_r),
									IDATA(freeLocal_u),
									IDATA(freeGlobal_u),
									IDATA(isDOFBoundary),
									sigma,
									DDATA(v),
									DDATA(n),
									DDATA(a),
									DDATA(grad_w),
									DDATA(dS),
									DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateInteriorElementBoundaryDiffusionAdjoint_sd(PyObject* self, 
							   PyObject* args)
{
  double sigma;
  PyObject *rowptr, *colind, *interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*u,*n,*a,*grad_w,*dS,*residual;
  if(!PyArg_ParseTuple(args,"OOOOOdOOOOOO",
		       &rowptr,
		       &colind,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &sigma,
		       &u,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
		       &residual))
    return NULL;
  updateInteriorElementBoundaryDiffusionAdjoint_sd(SHAPE(interiorElementBoundaries)[0],
						   SHAPE(grad_w)[1],
						   SHAPE(grad_w)[2],
						   SHAPE(grad_w)[3],
						   SHAPE(grad_w)[4],
						   IDATA(rowptr),
						   IDATA(colind),
						   IDATA(interiorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   sigma,
						   DDATA(u),
						   DDATA(n),
						   DDATA(a),
						   DDATA(grad_w),
						   DDATA(dS),
						   DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateExteriorElementBoundaryDiffusionAdjoint_sd(PyObject* self, 
							   PyObject* args)
{
  double sigma;
  PyObject *rowptr,*colind,*isDOFBoundary,*exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,*u,*ub,*n,*a,*grad_w,*dS,*residual;
  if(!PyArg_ParseTuple(args,"OOOOOOdOOOOOOO",
		       &rowptr,
		       &colind,
		       &isDOFBoundary,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &sigma,
		       &u,
		       &ub,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
		       &residual))
    return NULL;
  updateExteriorElementBoundaryDiffusionAdjoint_sd(SHAPE(exteriorElementBoundaries)[0],
						   SHAPE(grad_w)[1],
						   SHAPE(grad_w)[2],
						   SHAPE(grad_w)[3],
						   IDATA(rowptr),
						   IDATA(colind),
						   IDATA(isDOFBoundary),
						   IDATA(exteriorElementBoundaries),
						   IDATA(elementBoundaryElements),
						   IDATA(elementBoundaryLocalElementBoundaries),
						   sigma,
						   DDATA(u),
						   DDATA(ub),
						   DDATA(n),
						   DDATA(a),
						   DDATA(grad_w),
						   DDATA(dS),
						   DDATA(residual));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense_sd(PyObject* self, 
										   PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  double sigma;
  PyObject *rowptr,*colind,*interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *v,*n,*a,*grad_w,*dS,*jac;
  if(!PyArg_ParseTuple(args,"OOiiiiiOOOOOOOOOdOOOOOO",
                       &rowptr,
		       &colind,
		       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &sigma,
		       &v,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense_sd(SHAPE(interiorElementBoundaries)[0],
									   SHAPE(grad_w)[1],
									   SHAPE(grad_w)[2],
									   SHAPE(grad_w)[3],
									   SHAPE(v)[3],
									   SHAPE(n)[3],
									   IDATA(rowptr),
									   IDATA(colind),
									   offset_r,
									   stride_r,
									   offset_u,
									   stride_u,
									   nFreeVDOF_global,
									   IDATA(interiorElementBoundaries),
									   IDATA(elementBoundaryElements),
									   IDATA(elementBoundaryLocalElementBoundaries),
									   IDATA(nFreeDOF_element_r),
									   IDATA(nFreeDOF_element_u),
									   IDATA(freeLocal_r),
									   IDATA(freeGlobal_r),
									   IDATA(freeLocal_u),
									   IDATA(freeGlobal_u),
									   sigma,
									   DDATA(v),
									   DDATA(n),
									   DDATA(a),
									   DDATA(grad_w),
									   DDATA(dS),
									   DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense_sd(PyObject* self, 
										   PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  double sigma;
  PyObject *rowptr,*colind,*exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *isDOFBoundary,*v,*n,*a,*grad_w,*dS,*jac;
  if(!PyArg_ParseTuple(args,"OOiiiiiOOOOOOOOOOdOOOOOO",
                       &rowptr,
		       &colind,
		       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &isDOFBoundary,
		       &sigma,
		       &v,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense_sd(SHAPE(exteriorElementBoundaries)[0],
									   SHAPE(grad_w)[1],
									   SHAPE(grad_w)[2],
									   SHAPE(v)[2],
									   SHAPE(n)[2],
									   IDATA(rowptr),
									   IDATA(colind),
									   offset_r,
									   stride_r,
									   offset_u,
									   stride_u,
									   nFreeVDOF_global,
									   IDATA(exteriorElementBoundaries),
									   IDATA(elementBoundaryElements),
									   IDATA(elementBoundaryLocalElementBoundaries),
									   IDATA(nFreeDOF_element_r),
									   IDATA(nFreeDOF_element_u),
									   IDATA(freeLocal_r),
									   IDATA(freeGlobal_r),
									   IDATA(freeLocal_u),
									   IDATA(freeGlobal_u),
									   IDATA(isDOFBoundary),
									   sigma,
									   DDATA(v),
									   DDATA(n),
									   DDATA(a),
									   DDATA(grad_w),
									   DDATA(dS),
									   DDATA(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_CSR_sd(PyObject* self, 
										   PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  double sigma;
  PyObject *csrRowIndeces_ru,*csrColumnOffsets_eb_ru,*rowptr,*colind,*interiorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *v,*n,*a,*grad_w,*dS,*jac;
  if(!PyArg_ParseTuple(args,"OOiiiiiOOOOOOOOOOOdOOOOOO",
                       &rowptr,
		       &colind,
		       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &interiorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &csrRowIndeces_ru,
		       &csrColumnOffsets_eb_ru,
		       &sigma,
		       &v,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_CSR_sd(SHAPE(interiorElementBoundaries)[0],
									 SHAPE(grad_w)[1],
									 SHAPE(grad_w)[2],
									 SHAPE(grad_w)[3],
									 SHAPE(v)[3],
									 SHAPE(n)[3],
									 IDATA(rowptr),
									 IDATA(colind),
									 offset_r,
									 stride_r,
									 offset_u,
									 stride_u,
									 nFreeVDOF_global,
									 IDATA(interiorElementBoundaries),
									 IDATA(elementBoundaryElements),
									 IDATA(elementBoundaryLocalElementBoundaries),
									 IDATA(nFreeDOF_element_r),
									 IDATA(nFreeDOF_element_u),
									 IDATA(freeLocal_r),
									 IDATA(freeGlobal_r),
									 IDATA(freeLocal_u),
									 IDATA(freeGlobal_u),
									 IDATA(csrRowIndeces_ru),
									 IDATA(csrColumnOffsets_eb_ru),
									 sigma,
									 DDATA(v),
									 DDATA(n),
									 DDATA(a),
									 DDATA(grad_w),
									 DDATA(dS),
									 CSRVAL(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_CSR_sd(PyObject* self, 
										   PyObject* args)
{
  int offset_r,stride_r,offset_u,stride_u,nFreeVDOF_global;
  double sigma;
  PyObject *csrRowIndeces_ru,*csrColumnOffsets_eb_ru,*rowptr,*colind,*exteriorElementBoundaries,*elementBoundaryElements,*elementBoundaryLocalElementBoundaries,
    *nFreeDOF_element_r,*freeLocal_r,*freeGlobal_r,
    *nFreeDOF_element_u,*freeLocal_u,*freeGlobal_u,
    *isDOFBoundary,*v,*n,*a,*grad_w,*dS,*jac;
  if(!PyArg_ParseTuple(args,"OOiiiiiOOOOOOOOOOOOdOOOOOO",
                       &rowptr,
		       &colind,
		       &offset_r,
                       &stride_r,
                       &offset_u,
                       &stride_u,
		       &nFreeVDOF_global,
                       &exteriorElementBoundaries,
                       &elementBoundaryElements,
                       &elementBoundaryLocalElementBoundaries,
                       &nFreeDOF_element_r,
                       &freeLocal_r,
                       &freeGlobal_r,
                       &nFreeDOF_element_u,
                       &freeLocal_u,
                       &freeGlobal_u,
		       &csrRowIndeces_ru,
		       &csrColumnOffsets_eb_ru,
		       &isDOFBoundary,
		       &sigma,
		       &v,
		       &n,
		       &a,
		       &grad_w,
		       &dS,
                       &jac))
    return NULL;
  updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_CSR_sd(SHAPE(exteriorElementBoundaries)[0],
									 SHAPE(grad_w)[1],
									 SHAPE(grad_w)[2],
									 SHAPE(v)[2],
									 SHAPE(n)[2],
									 IDATA(rowptr),
									 IDATA(colind),
									 offset_r,
									 stride_r,
									 offset_u,
									 stride_u,
									 nFreeVDOF_global,
									 IDATA(exteriorElementBoundaries),
									 IDATA(elementBoundaryElements),
									 IDATA(elementBoundaryLocalElementBoundaries),
									 IDATA(nFreeDOF_element_r),
									 IDATA(nFreeDOF_element_u),
									 IDATA(freeLocal_r),
									 IDATA(freeGlobal_r),
									 IDATA(freeLocal_u),
									 IDATA(freeGlobal_u),
									 IDATA(csrRowIndeces_ru),
									 IDATA(csrColumnOffsets_eb_ru),
									 IDATA(isDOFBoundary),
									 sigma,
									 DDATA(v),
									 DDATA(n),
									 DDATA(a),
									 DDATA(grad_w),
									 DDATA(dS),
									 CSRVAL(jac));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdate_f_movingDomain(PyObject* self, 
				   PyObject* args)
{
  PyObject *xt,*m,*f;
  if(!PyArg_ParseTuple(args,"OOO",
                       &xt,
		       &m,
		       &f))
    return NULL;
  if (ND(m) == 2)//we can use q for ebqe too
    update_f_movingDomain_q(SHAPE(f)[0],
			    SHAPE(f)[1],
			    SHAPE(f)[2],
			    DDATA(xt),
			    DDATA(m),
			    DDATA(f));
  else if (ND(m) == 3)
    update_f_movingDomain_ebq(SHAPE(f)[0],
			      SHAPE(f)[1],
			      SHAPE(f)[2],
			      SHAPE(f)[3],
			      DDATA(xt),
			      DDATA(m),
			      DDATA(f));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdate_f_movingDomain_constantMass(PyObject* self, 
						PyObject* args)
{
  PyObject *xt,*f;
  if(!PyArg_ParseTuple(args,"OO",
                       &xt,
		       &f))
    return NULL;
  if (ND(f) == 3)//we can use q for ebqe too
    update_f_movingDomain_constantMass_q(SHAPE(f)[0],
					 SHAPE(f)[1],
					 SHAPE(f)[2],
					 DDATA(xt),
					 DDATA(f));
  else if (ND(f) == 4)
    update_f_movingDomain_constantMass_ebq(SHAPE(f)[0],
					   SHAPE(f)[1],
					   SHAPE(f)[2],
					   SHAPE(f)[3],
					   DDATA(xt),
					   DDATA(f));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
cfemIntegralsUpdate_df_movingDomain(PyObject* self, 
				   PyObject* args)
{
  PyObject *xt,*dm,*df;
  if(!PyArg_ParseTuple(args,"OOO",
                       &xt,
		       &dm,
		       &df))
    return NULL;
  //we can use the update_f functions for this
  if (ND(dm) == 2) //we can use _q for _ebqe too
    update_f_movingDomain_q(SHAPE(df)[0],
			    SHAPE(df)[1],
			    SHAPE(df)[2],
			    DDATA(xt),
			    DDATA(dm),
			    DDATA(df));
  else if (ND(dm) == 3)
    update_f_movingDomain_ebq(SHAPE(df)[0],
			      SHAPE(df)[1],
			      SHAPE(df)[2],
			      SHAPE(df)[3],
			      DDATA(xt),
			      DDATA(dm),
			      DDATA(df));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateStress_weak(PyObject* self, 
						PyObject* args)
{
  PyObject *sigma,
    *grad_w_dV,
    *weak_residual_x,
    *weak_residual_y,
    *weak_residual_z;
  if(!PyArg_ParseTuple(args,"OOOOO",
		       &sigma,
		       &grad_w_dV,
		       &weak_residual_x,
		       &weak_residual_y,
		       &weak_residual_z))
    return NULL;
  updateStress_weak(SHAPE(grad_w_dV)[0],
		    SHAPE(grad_w_dV)[1],
		    SHAPE(grad_w_dV)[2],
		    SHAPE(grad_w_dV)[3],
		    DDATA(sigma),
		    DDATA(grad_w_dV),
		    DDATA(weak_residual_x),
		    DDATA(weak_residual_y),
		    DDATA(weak_residual_z));
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* cfemIntegralsUpdateStressJacobian_weak(PyObject* self, 
							PyObject* args)
{
  PyObject *dsigma_xx,
    *dsigma_xy,
    *dsigma_xz,
    *dsigma_yx,
    *dsigma_yy,
    *dsigma_yz,
    *dsigma_zx,
    *dsigma_zy,
    *dsigma_zz,
    *grad_v,
    *grad_w_dV,
    *jacobian_weak_residual_xx,
    *jacobian_weak_residual_xy,
    *jacobian_weak_residual_xz,
    *jacobian_weak_residual_yx,
    *jacobian_weak_residual_yy,
    *jacobian_weak_residual_yz,
    *jacobian_weak_residual_zx,
    *jacobian_weak_residual_zy,
    *jacobian_weak_residual_zz;
  if(!PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOO",
		       &dsigma_xx,
		       &dsigma_xy,
		       &dsigma_xz,
		       &dsigma_yx,
		       &dsigma_yy,
		       &dsigma_yz,
		       &dsigma_zx,
		       &dsigma_zy,
		       &dsigma_zz,
		       &grad_v,
		       &grad_w_dV,
		       &jacobian_weak_residual_xx,
		       &jacobian_weak_residual_xy,
		       &jacobian_weak_residual_xz,
		       &jacobian_weak_residual_yx,
		       &jacobian_weak_residual_yy,
		       &jacobian_weak_residual_yz,
		       &jacobian_weak_residual_zx,
		       &jacobian_weak_residual_zy,
		       &jacobian_weak_residual_zz))
    return NULL;
  updateStressJacobian_weak(SHAPE(grad_w_dV)[0],
			    SHAPE(grad_w_dV)[1],
			    SHAPE(grad_w_dV)[2],
			    SHAPE(grad_w_dV)[2],
			    SHAPE(grad_w_dV)[3],
			    DDATA(dsigma_xx),
			    DDATA(dsigma_xy),
			    DDATA(dsigma_xz),
			    DDATA(dsigma_yx),
			    DDATA(dsigma_yy),
			    DDATA(dsigma_yz),
			    DDATA(dsigma_zx),
			    DDATA(dsigma_zy),
			    DDATA(dsigma_zz),
			    DDATA(grad_v),
			    DDATA(grad_w_dV),
			    DDATA(jacobian_weak_residual_xx),
			    DDATA(jacobian_weak_residual_xy),
			    DDATA(jacobian_weak_residual_xz),
			    DDATA(jacobian_weak_residual_yx),
			    DDATA(jacobian_weak_residual_yy),
			    DDATA(jacobian_weak_residual_yz),
			    DDATA(jacobian_weak_residual_zx),
			    DDATA(jacobian_weak_residual_zy),
			    DDATA(jacobian_weak_residual_zz));
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* cfemIntegralsProjectFromNodalInterpolationConditions(PyObject* self, 
								       PyObject* args)
{
  int dim_dof;
  PyObject *l2g,
    *functional_map_element,
    *interpolationValues,
    *dofs;

  if (!PyArg_ParseTuple(args,"iOOOO",
			&dim_dof,
			&l2g,
			&functional_map_element,
			&interpolationValues,
			&dofs))
    return NULL;
  
  projectFromNodalInterpolationConditions(SHAPE(l2g)[0],
					  SHAPE(l2g)[1],
					  dim_dof,
					  IDATA(l2g),
					  IDATA(functional_map_element),
					  DDATA(interpolationValues),
					  DDATA(dofs));

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef cfemIntegralsMethods[] = {
 { "calculateVelocityProjectionMatrixLDG",
    cfemIntegralsCalculateVelocityProjectionMatrixLDG,
    METH_VARARGS, 
    "Calculate the projection matrix for obtaining the diffusive velocity degrees of freedom"},
  { "calculateVelocityQuadrature_MixedForm_Jacobian",
    cfemIntegralsCalculateVelocityQuadrature_MixedForm_Jacobian,
    METH_VARARGS, 
    "Calculate the Jacobian of the diffusive velocity at the quadrature points for mixed form diffusion"},
  { "calculateVelocityQuadrature_MixedForm",
    cfemIntegralsCalculateVelocityQuadrature_MixedForm,
    METH_VARARGS, 
    "Calculate the diffusive velocity at the quadrature points for mixed form diffusion"},
  { "calculateVelocityQuadrature_MixedForm2_Jacobian",
    cfemIntegralsCalculateVelocityQuadrature_MixedForm2_Jacobian,
    METH_VARARGS, 
    "Calculate the Jacobian of the diffusive velocity at the quadrature points for mixed form diffusion"},
  { "calculateVelocityQuadrature_MixedForm2_Jacobian_sd",
    cfemIntegralsCalculateVelocityQuadrature_MixedForm2_Jacobian_sd,
    METH_VARARGS, 
    "Calculate the Jacobian of the diffusive velocity at the quadrature points for mixed form diffusion sparse rep"},
  { "calculateVelocityQuadrature_MixedForm2",
    cfemIntegralsCalculateVelocityQuadrature_MixedForm2,
    METH_VARARGS, 
    "Calculate the diffusive velocity at the quadrature points for mixed form diffusion"},
  { "calculateVelocityQuadrature_MixedForm2_sd",
    cfemIntegralsCalculateVelocityQuadrature_MixedForm2_sd,
    METH_VARARGS, 
    "Calculate the diffusive velocity at the quadrature points for mixed form diffusion sparse rep"},
  { "calculateVelocityQuadrature_MixedForm2_vdof_sd",
    cgwvdNumericalFluxCalculateVelocityQuadrature_MixedForm2_vdof_sd,
    METH_VARARGS, 
    "Calculate the diffusive velocity at the quadrature points for mixed form diffusion sparse rep and saving the velocity degrees of freedom"},
  { "updatePotential_MixedForm_weakJacobian",
    cfemIntegralsUpdatePotential_MixedForm_weakJacobian,
    METH_VARARGS, 
    "Update the element Jacobian integral for the mixed form of diffusion"},
  { "updatePotential_MixedForm_weak",
    cfemIntegralsUpdatePotential_MixedForm_weak,
    METH_VARARGS, 
    "Update the element integral for the mixed form of diffusion"},
 { "updatePotential_MixedForm_weak_gwvd",
    cfemIntegralsUpdatePotential_MixedForm_weak_gwvd,
    METH_VARARGS, 
    "Update the element integral for the mixed form of diffusion for variable density flow LDG Method"},
  { "updateExteriorElementBoundary_MixedForm_weakJacobian",
    cfemIntegralsUpdateExteriorElementBoundary_MixedForm_weakJacobian,
    METH_VARARGS, 
    "Update the Jacobian of the boundary term for the mixed form diffusion on exterior element boundaries"}, 
  { "updateExteriorElementBoundary_MixedForm_weak",
    cfemIntegralsUpdateExteriorElementBoundary_MixedForm_weak,
    METH_VARARGS, 
    "Update the boundary term for the mixed form diffusion on exterior element boundaries"}, 
  { "updateInteriorElementBoundary_MixedForm_weakJacobian",
    cfemIntegralsUpdateInteriorElementBoundary_MixedForm_weakJacobian,
    METH_VARARGS, 
    "Update the Jacobian of the boundary term for the mixed form diffusion on interior element boundaries"}, 
  { "updateInteriorElementBoundary_MixedForm_weak",
    cfemIntegralsUpdateInteriorElementBoundary_MixedForm_weak,
    METH_VARARGS, 
    "Update the boundary term for mixed form diffusion on interior element boundaries"}, 
  { "calculateExteriorNumericalTrace_Potential",
    cfemIntegralsCalculateExteriorNumericalTrace_Potential,
    METH_VARARGS, 
    "Calcualte the trace of the potential on interior element boundaries"}, 
  { "calculateInteriorNumericalTrace_Potential",
    cfemIntegralsCalculateInteriorNumericalTrace_Potential,
    METH_VARARGS, 
    "Calculate the trace of the potential on interior element boundaries"}, 
  { "updateMass_weak", 
    cfemIntegralsUpdateMass_weak, 
    METH_VARARGS, 
    "Update the weak residual with the weak mass integral"}, 
  { "updateMassJacobian_weak", 
    cfemIntegralsUpdateMassJacobian_weak, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak mass integral"}, 
  { "updateMassJacobian_weak_lowmem", 
    cfemIntegralsUpdateMassJacobian_weak_lowmem, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak mass integral"}, 
  { "updateMass_strong", 
    cfemIntegralsUpdateMass_strong, 
    METH_VARARGS, 
    "Update the strong residual with the strong mass term"}, 
  { "updateMassJacobian_strong", 
    cfemIntegralsUpdateMassJacobian_strong, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong mass term"}, 
  { "updateMass_adjoint", 
    cfemIntegralsUpdateMass_adjoint, 
    METH_VARARGS, 
    "Update the adjoint applied to the test  functions with the mass term"}, 
  { "updateAdvection_weak", 
    cfemIntegralsUpdateAdvection_weak, 
    METH_VARARGS, 
    "Update the weak residual with the weak advection integral"}, 
  { "updateAdvectionJacobian_weak", 
    cfemIntegralsUpdateAdvectionJacobian_weak, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak advection integral"}, 
  { "updateAdvectionJacobian_weak_lowmem", 
    cfemIntegralsUpdateAdvectionJacobian_weak_lowmem, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak advection integral"}, 
  { "updateAdvection_strong", 
    cfemIntegralsUpdateAdvection_strong, 
    METH_VARARGS, 
    "Update the strong residual with the strong advection term"}, 
  { "updateAdvectionJacobian_strong", 
    cfemIntegralsUpdateAdvectionJacobian_strong, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong advection term"}, 
  { "updateAdvection_adjoint", 
    cfemIntegralsUpdateAdvection_adjoint, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the advection term"}, 
  { "updateHamiltonian_weak", 
    cfemIntegralsUpdateHamiltonian_weak, 
    METH_VARARGS, 
    "Update the weak residual with the weak Hamiltonian integral"}, 
  { "updateHamiltonianJacobian_weak", 
    cfemIntegralsUpdateHamiltonianJacobian_weak, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak Hamiltonian integral"}, 
  { "updateHamiltonianJacobian_weak_lowmem", 
    cfemIntegralsUpdateHamiltonianJacobian_weak_lowmem, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak Hamiltonian integral"}, 
  { "updateHamiltonian_strong", 
    cfemIntegralsUpdateHamiltonian_strong, 
    METH_VARARGS, 
    "Update the strong residual with the strong Hamiltonian term"}, 
  { "updateHamiltonianJacobian_strong", 
    cfemIntegralsUpdateHamiltonianJacobian_strong, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong Hamiltonian term"}, 
  { "updateHamiltonian_adjoint", 
    cfemIntegralsUpdateHamiltonian_adjoint, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the Hamiltonian term"}, 
  { "updateDiffusion_weak", 
    cfemIntegralsUpdateDiffusion_weak, 
    METH_VARARGS, 
    "Update the weak residual with the weak diffusion integral"}, 
  { "updateDiffusion_weak_lowmem", 
    cfemIntegralsUpdateDiffusion_weak_lowmem, 
    METH_VARARGS, 
    "Update the weak residual with the weak diffusion integral"}, 
  { "updateDiffusion_weak_sd", 
    cfemIntegralsUpdateDiffusion_weak_sd, 
    METH_VARARGS, 
    "Update the weak residual with the weak diffusion integral (sparse diffusion tensor)"}, 
  { "updateDiffusionJacobian_weak", 
    cfemIntegralsUpdateDiffusionJacobian_weak, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak diffusion integral"}, 
  { "updateDiffusionJacobian_weak_lowmem", 
    cfemIntegralsUpdateDiffusionJacobian_weak_lowmem, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak diffusion integral"}, 
  { "updateDiffusionJacobian_weak_sd", 
    cfemIntegralsUpdateDiffusionJacobian_weak_sd, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak diffusion integral (sparse diffusion tensor)"}, 
  { "updateDiffusion_MixedForm_weak", 
    cfemIntegralsUpdateDiffusion_MixedForm_weak, 
    METH_VARARGS, 
    "Update the weak residual with the weak diffusion integral"}, 
  { "updateDiffusion_MixedForm_weak_sd", 
    cfemIntegralsUpdateDiffusion_MixedForm_weak_sd, 
    METH_VARARGS, 
    "Update the weak residual with the weak diffusion integral (sparse diffusion tensor)"}, 
  { "updateDiffusionJacobian_MixedForm_weak", 
    cfemIntegralsUpdateDiffusionJacobian_MixedForm_weak, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak diffusion integral"}, 
  { "updateDiffusionJacobian_MixedForm_weak_sd", 
    cfemIntegralsUpdateDiffusionJacobian_MixedForm_weak_sd, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak diffusion integral (sparse diffusion tensor)"}, 
  { "updateDiffusion_strong", 
    cfemIntegralsUpdateDiffusion_strong, 
    METH_VARARGS, 
    "Update the strong residual with the strong diffusion term"}, 
  { "updateDiffusion_strong_sd", 
    cfemIntegralsUpdateDiffusion_strong_sd, 
    METH_VARARGS, 
    "Update the strong residual with the strong diffusion term (sparse diffusion tensor)"}, 
  { "updateDiffusionJacobian_strong", 
    cfemIntegralsUpdateDiffusionJacobian_strong, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong diffusion term"}, 
  { "updateDiffusionJacobian_strong_sd", 
    cfemIntegralsUpdateDiffusionJacobian_strong_sd, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong diffusion term (sparse diffusion tensor)"}, 
  { "updateDiffusion2_strong", 
    cfemIntegralsUpdateDiffusion2_strong, 
    METH_VARARGS, 
    "Update the strong residual with the strong diffusion term"}, 
  { "updateDiffusion2_strong_sd", 
    cfemIntegralsUpdateDiffusion2_strong_sd, 
    METH_VARARGS, 
    "Update the strong residual with the strong diffusion term (sparse diffusion tensor)"}, 
  { "updateDiffusionJacobian2_strong", 
    cfemIntegralsUpdateDiffusionJacobian2_strong, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong diffusion term"}, 
  { "updateDiffusionJacobian2_strong_sd", 
    cfemIntegralsUpdateDiffusionJacobian2_strong_sd, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong diffusion term (sparse diffusion tensor)"}, 
  { "updateDiffusion_adjoint", 
    cfemIntegralsUpdateDiffusion_adjoint, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the diffusion term"}, 
  { "updateDiffusion_adjoint_sd", 
    cfemIntegralsUpdateDiffusion_adjoint_sd, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the diffusion term sparse diffusion tensor"}, 
  { "updateDiffusion2_adjoint", 
    cfemIntegralsUpdateDiffusion2_adjoint, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the diffusion term"}, 
  { "updateDiffusion2_adjoint_sd", 
    cfemIntegralsUpdateDiffusion2_adjoint_sd, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the diffusion term sparse diffusion tensor"}, 
  { "updateReaction_weak", 
    cfemIntegralsUpdateReaction_weak, 
    METH_VARARGS, 
    "Update the weak residual with the weak reaction integral"}, 
  { "updateReactionJacobian_weak", 
    cfemIntegralsUpdateReactionJacobian_weak, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak reaction integral"}, 
  { "updateReactionJacobian_weak_lowmem", 
    cfemIntegralsUpdateReactionJacobian_weak_lowmem, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak reaction integral"}, 
  { "updateReaction_strong", 
    cfemIntegralsUpdateReaction_strong, 
    METH_VARARGS, 
    "Update the strong residual with the strong reaction term"}, 
  { "updateReactionJacobian_strong", 
    cfemIntegralsUpdateReactionJacobian_strong, 
    METH_VARARGS, 
    "Update the Jacobian of the strong residual with the Jacobian of the strong reaction term"}, 
  { "updateReaction_adjoint", 
    cfemIntegralsUpdateReaction_adjoint, 
    METH_VARARGS, 
    "Update the adjoint applied to the test functions with the reaction term"}, 
  { "updateSubgridError", 
    cfemIntegralsUpdateSubgridError, 
    METH_VARARGS, 
    "Update the weak residual with the weak subgrid error integral"}, 
  { "updateSubgridErrorJacobian", 
    cfemIntegralsUpdateSubgridErrorJacobian, 
    METH_VARARGS, 
    "Update Jacobian of the weak residual with the Jacobian of the weak subgrid error integral"}, 
  { "updateNumericalDiffusion", 
    cfemIntegralsUpdateNumericalDiffusion, 
    METH_VARARGS, 
    "Update the weak residual with the weak numerical diffusion integral"}, 
  { "updateNumericalDiffusion_lowmem", 
    cfemIntegralsUpdateNumericalDiffusion_lowmem, 
    METH_VARARGS, 
    "Update the weak residual with the weak numerical diffusion integral"}, 
  { "updateNumericalDiffusionJacobian", 
    cfemIntegralsUpdateNumericalDiffusionJacobian, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak numerical diffusion integral"}, 
  { "updateNumericalDiffusionJacobian_lowmem", 
    cfemIntegralsUpdateNumericalDiffusionJacobian_lowmem, 
    METH_VARARGS, 
    "Update the Jacobian of the weak residual with the Jacobian of the weak numerical diffusion integral"}, 
  { "calculateScalarScalarProduct", 
    cfemIntegralsCalculateScalarScalarProduct, 
    METH_VARARGS, 
    "Scalar multiplication of a scalar"},
  { "calculateVectorScalarProduct", 
    cfemIntegralsCalculateVectorScalarProduct, 
    METH_VARARGS, 
    "Scalar multiplication of a vector"},
  { "calculateTensorScalarProduct", 
    cfemIntegralsCalculateTensorScalarProduct, 
    METH_VARARGS, 
    "Scalar multiplication of a tensor"},
  { "calculateWeightedShape", 
    cfemIntegralsCalculateWeightedShape,
    METH_VARARGS, 
    "The quadrature weighted shape  functions"},
  { "calculateWeightedShapeGradients", 
    cfemIntegralsCalculateWeightedShapeGradients,
    METH_VARARGS, 
    "The quadrature weighted shape  function gradients"},
  { "calculateWeightedPiolaShapeGradients", 
    cfemIntegralsCalculateWeightedPiolaShapeGradients,
    METH_VARARGS, 
    "The quadrature weighted shape function gradients using Piola Mapping"},
  { "calculateWeightedShapeHessians", 
    cfemIntegralsCalculateWeightedShapeHessians,
    METH_VARARGS, 
    "The quadrature weighted shape  function gradients"},
  { "calculateShape_X_weightedShape", 
    cfemIntegralsCalculateShape_X_weightedShape,
    METH_VARARGS, 
    "The tensor product of  shape  functions"},
  { "calculateShape_X_weightedGradShape", 
    cfemIntegralsCalculateShape_X_weightedGradShape,
    METH_VARARGS, 
    "The tensor product of  shape and shape  gradient functions"},
  { "calculateGradShape_X_weightedShape", 
    cfemIntegralsCalculateGradShape_X_weightedShape,
    METH_VARARGS, 
    "The tensor product of  shape gradient and shape functions"},
  { "calculateGradShape_X_weightedGradShape", 
    cfemIntegralsCalculateGradShape_X_weightedGradShape,
    METH_VARARGS, 
    "The tensor product of  shape gradient and shape gradient functions"},
  { "calculateWeightedShapeTrace", 
    cfemIntegralsCalculateWeightedShapeTrace,
    METH_VARARGS, 
    "The quadrature weighted shape  functions"},
  { "calculateWeightedPiolaShapeTrace", 
    cfemIntegralsCalculateWeightedPiolaShapeTrace,
    METH_VARARGS, 
    "The quadrature weighted shape functions for Piola mapping"},
  { "calculateShape_X_weightedShapeTrace", 
    cfemIntegralsCalculateShape_X_weightedShapeTrace,
    METH_VARARGS, 
    "The tensor product of  shape  functions"},
  { "calculateGradShape_X_weightedShapeTrace", 
    cfemIntegralsCalculateGradShape_X_weightedShapeTrace,
    METH_VARARGS, 
    "The tensor product of  shape gradient and shape functions"},
  { "calculateWeightedShapeGlobalExteriorTrace", 
    cfemIntegralsCalculateWeightedShapeGlobalExteriorTrace,
    METH_VARARGS, 
    "The quadrature weighted shape  functions"},
  { "calculateShape_X_weightedShapeGlobalExteriorTrace", 
    cfemIntegralsCalculateShape_X_weightedShapeGlobalExteriorTrace,
    METH_VARARGS, 
    "The tensor product of  shape  functions"},
  { "calculateGradShape_X_weightedShapeGlobalExteriorTrace", 
    cfemIntegralsCalculateGradShape_X_weightedShapeGlobalExteriorTrace,
    METH_VARARGS, 
    "The tensor product of  shape gradient and shape functions"},
  { "calculateFiniteElementFunctionValues", 
    cfemIntegralsCalculateFiniteElementFunctionValues,
    METH_VARARGS, 
    "Calculate the values of the finite element function at quadrature points"},
  { "calculateFiniteElementFunctionGradientValues", 
    cfemIntegralsCalculateFiniteElementFunctionGradientValues,
    METH_VARARGS, 
    "Calculate the values of the gradients of the finite element function at quadrature points"},
  { "calculateFiniteElementFunctionHessianValues", 
    cfemIntegralsCalculateFiniteElementFunctionHessianValues,
    METH_VARARGS, 
    "Calculate the values of the gradients of the finite element function at quadrature points"},
  { "calculateFiniteElementFunctionValuesTrace", 
    cfemIntegralsCalculateFiniteElementFunctionValuesTrace,
    METH_VARARGS, 
    "Calculate the values of the finite element function at quadrature points on the element boundary"},
  { "calculateFiniteElementFunctionGradientValuesTrace", 
    cfemIntegralsCalculateFiniteElementFunctionGradientValuesTrace,
    METH_VARARGS, 
    "Calculate the values of the gradients of the finite element function at quadrature points on the element boundary"},
  { "calculateFiniteElementFunctionValuesGlobalExteriorTrace", 
    cfemIntegralsCalculateFiniteElementFunctionValuesGlobalExteriorTrace,
    METH_VARARGS, 
    "Calculate the values of the finite element function at quadrature points on the global boundary"},
  { "calculateFiniteElementFunctionGradientValuesGlobalExteriorTrace", 
    cfemIntegralsCalculateFiniteElementFunctionGradientValuesGlobalExteriorTrace,
    METH_VARARGS, 
    "Calculate the values of the gradients of the finite element function at quadrature points on the global boundary"},
  { "calculateFiniteElementFunctionGradientTensorValues", 
    cfemIntegralsCalculateFiniteElementFunctionGradientTensorValues, 
    METH_VARARGS, 
    "Calculate tensor product  of finite element function gradients with test function gradients"},
  { "updateGlobalResidualFromElementResidual", 
    cfemIntegralsUpdateGlobalResidualFromElementResidual,
    METH_VARARGS, 
    "load  the global residual"},
  { "updateGlobalJacobianFromElementJacobian_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromElementJacobian_CSR,
    METH_VARARGS, 
    "load  the global jacobian"},
  { "updateGlobalJacobianFromElementJacobian_eb_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromElementJacobian_eb_CSR,
    METH_VARARGS, 
    "load  the global jacobian"},
  { "updateGlobalJacobianFromElementJacobian_dense", 
    cfemIntegralsUpdateGlobalJacobianFromElementJacobian_dense,
    METH_VARARGS, 
    "load  the global jacobian"},
  { "updateGlobalJacobianFromElementJacobian_eb_dense", 
    cfemIntegralsUpdateGlobalJacobianFromElementJacobian_eb_dense,
    METH_VARARGS, 
    "load  the global jacobian"},
  { "calculateFlowVelocity", 
    cfemIntegralsCalculateFlowVelocity,
    METH_VARARGS, 
    "calculate the total flow velocity"},
  { "updateAddJacobian_CSR", 
    cfemIntegralsUpdateAddJacobian_CSR,
    METH_VARARGS, 
    "add a jacobian entry"},
  { "zeroJacobian_CSR", 
    cfemIntegralsZeroJacobian_CSR,
    METH_VARARGS, 
    "zero a jacobian entry"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the external boundary"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the external boundary"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_CSR", 
    cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_CSR,
    METH_VARARGS, 
    "set the  jacobian entries due  to the Hamilton Jacobi numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense", 
    cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_dense", 
    cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_eb_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the internal boundaries"},
  { "updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense", 
    cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the external boundary"},
  { "updateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_dense", 
    cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryFluxJacobian_eb_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the numerical flux on the external boundary"},
  { "updateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_dense", 
    cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryFluxJacobian_2sided_dense,
    METH_VARARGS, 
    "set the  jacobian entries due  to the Hamilton Jacobi numerical flux on the internal boundaries"},
  {"updateInteriorElementBoundaryFlux",
   cfemIntegralsUpdateInteriorElementBoundaryFlux,
   METH_VARARGS,
   "update the residual with the element boundary fluxes on the interior element boundaries"},
  {"updateExteriorElementBoundaryFlux",
   cfemIntegralsUpdateExteriorElementBoundaryFlux,
   METH_VARARGS,
   "update the residual with the element boundary fluxes on the exterior element boundaries"},
  {"updateExteriorElementBoundaryStressFlux",
   cfemIntegralsUpdateExteriorElementBoundaryStressFlux,
   METH_VARARGS,
   "update the residual with the element boundary fluxes on the exterior element boundaries"},
  {"updateInteriorTwoSidedElementBoundaryFlux",
   cfemIntegralsUpdateTwoSidedInteriorElementBoundaryFlux,
   METH_VARARGS,
   "update the residual with the two-sided element boundary fluxes on the interior element boundaries"},
  {"calculateInteriorElementBoundaryVelocities",
   cfemIntegralsCalculateInteriorElementBoundaryVelocities,
   METH_VARARGS,
   "calcualte the averages and jumps of mass and velocity on the interior element boundaries"},
  {"calculateExteriorElementBoundaryVelocities",
   cfemIntegralsCalculateExteriorElementBoundaryVelocities,
   METH_VARARGS,
   "calcualte the averages and jumps of mass and velocity on the interior element boundaries"},
  {"writeDOF",
   cfemIntegralsWriteDOF,
   METH_VARARGS,
   "write a component of the DOF to a file using the given format"},
  {"writeDOF_ZEROS",
   cfemIntegralsWriteDOF_ZEROS,
   METH_VARARGS,
   "write a component of the DOF to a file using the given format"},
  {"parametricMaps_getPermutations",
   cfemIntegralsParametricMaps_getPermutations,
   METH_VARARGS,
   "need to add doc"},
 {"getPermutationsGlobal",
  cfemIntegralsGetPermutationsGlobal,
  METH_VARARGS,
  "need to add doc"},
/*   {"parametricMaps_getPermutationsGlobalExterior", */
/*    cfemIntegralsParametricMaps_getPermutationsGlobalExterior, */
/*    METH_VARARGS, */
/*    "need to add doc"}, */
  {"parametricMaps_getValues",
   cfemIntegralsParametricMaps_getValues,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getValuesTrace",
   cfemIntegralsParametricMaps_getValuesTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getValuesGlobalExteriorTrace",
   cfemIntegralsParametricMaps_getValuesGlobalExteriorTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getInverseValues",
   cfemIntegralsParametricMaps_getInverseValues,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getInverseValuesTrace",
   cfemIntegralsParametricMaps_getInverseValuesTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getInverseValuesGlobalExteriorTrace",
   cfemIntegralsParametricMaps_getInverseValuesGlobalExteriorTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getJacobianValues",
   cfemIntegralsParametricMaps_getJacobianValues,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getJacobianValuesTrace",
   cfemIntegralsParametricMaps_getJacobianValuesTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getJacobianValuesGlobalExteriorTrace",
   cfemIntegralsParametricMaps_getJacobianValuesGlobalExteriorTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getJacobianValuesTrace_movingDomain",
   cfemIntegralsParametricMaps_getJacobianValuesTrace_movingDomain,
   METH_VARARGS,
   "need to add doc"},
  {"parametricMaps_getJacobianValuesGlobalExteriorTrace_movingDomain",
   cfemIntegralsParametricMaps_getJacobianValuesGlobalExteriorTrace_movingDomain,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getValues",
   cfemIntegralsParametricFiniteElementSpace_getValues,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getValuesTrace",
   cfemIntegralsParametricFiniteElementSpace_getValuesTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getValuesGlobalExteriorTrace",
   cfemIntegralsParametricFiniteElementSpace_getValuesGlobalExteriorTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getGradientValues",
   cfemIntegralsParametricFiniteElementSpace_getGradientValues,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getHessianValues",
   cfemIntegralsParametricFiniteElementSpace_getHessianValues,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getGradientValuesTrace",
   cfemIntegralsParametricFiniteElementSpace_getGradientValuesTrace,
   METH_VARARGS,
   "need to add doc"},
  {"parametricFiniteElementSpace_getGradientValuesGlobalExteriorTrace",
   cfemIntegralsParametricFiniteElementSpace_getGradientValuesGlobalExteriorTrace,
   METH_VARARGS,
   "need to add doc"},
  {"copyLeftElementBoundaryInfo",
   cfemIntegralsCopyLeftElementBoundaryInfo,
   METH_VARARGS,
   "need to add doc"},
  {"copyGlobalElementBoundaryInfo",
   cfemIntegralsCopyGlobalElementBoundaryInfo,
   METH_VARARGS,
   "need to add doc"},
  {"copyLeftElementBoundaryInfo_movingDomain",
   cfemIntegralsCopyLeftElementBoundaryInfo_movingDomain,
   METH_VARARGS,
   "need to add doc"},
  { "calculateDimensionlessNumbersADR", 
    cfemIntegralsCalculateDimensionlessNumbersADR,
    METH_VARARGS, 
    "Default manual calculation of CFL and Peclet numbers"},
  { "calculateDimensionlessNumbersADR_sd", 
    cfemIntegralsCalculateDimensionlessNumbersADR_sd,
    METH_VARARGS, 
    "Default manual calculation of CFL and Peclet numbers (sparse diffusion tensor)"},
  { "calculateCFLADR", 
    cfemIntegralsCalculateCFLADR,
    METH_VARARGS, 
    "Default manual calculation of CFL"},
  { "calculateCFLADR2speeds", 
    cfemIntegralsCalculateCFLADR2speeds,
    METH_VARARGS, 
    "Default manual calculation of CFL if have two speeds defined (not very realistic)"},
  { "calculateIntegrationWeights", 
    cfemIntegralsCalculateIntegrationWeights,
    METH_VARARGS, 
    "calculate integration weights in physical space"},
  { "calculateElementBoundaryIntegrationWeights", 
    cfemIntegralsCalculateElementBoundaryIntegrationWeights,
    METH_VARARGS, 
    "calculate integration weights for element boundaries in physical space"},
  { "updateInteriorElementBoundaryAdvectiveVelocity", 
    cfemIntegralsUpdateInteriorElementBoundaryAdvectiveVelocity,
    METH_VARARGS, 
    "add doc"},
  { "updateExteriorElementBoundaryAdvectiveVelocity", 
    cfemIntegralsUpdateExteriorElementBoundaryAdvectiveVelocity,
    METH_VARARGS, 
    "add doc"},
  { "updateInteriorElementBoundaryDiffusiveVelocity", 
    cfemIntegralsUpdateInteriorElementBoundaryDiffusiveVelocity,
    METH_VARARGS, 
    "add doc"},
  { "updateInteriorElementBoundaryDiffusiveVelocity_sd", 
    cfemIntegralsUpdateInteriorElementBoundaryDiffusiveVelocity_sd,
    METH_VARARGS, 
    "add doc"},
  { "updateExteriorElementBoundaryDiffusiveVelocity", 
    cfemIntegralsUpdateExteriorElementBoundaryDiffusiveVelocity,
    METH_VARARGS, 
    "add doc"},
  { "updateExteriorElementBoundaryDiffusiveVelocity_sd", 
    cfemIntegralsUpdateExteriorElementBoundaryDiffusiveVelocity_sd,
    METH_VARARGS, 
    "add doc"},
  { "calculateInteriorElementBoundaryAverageVelocity", 
    cfemIntegralsCalculateInteriorElementBoundaryAverageVelocity,
    METH_VARARGS, 
    "add doc"},
  { "updateInteriorElementBoundaryShockCapturingVelocity", 
    cfemIntegralsUpdateInteriorElementBoundaryShockCapturingVelocity,
    METH_VARARGS, 
    "mwf fooling around with post processing when have shock capturing too"},
  { "updateExteriorElementBoundaryShockCapturingVelocity", 
    cfemIntegralsUpdateExteriorElementBoundaryShockCapturingVelocity,
    METH_VARARGS, 
    "mwf fooling around with post processing when have shock capturing too"},
  { "calculateExteriorElementBoundaryAverageVelocity", 
    cfemIntegralsCalculateExteriorElementBoundaryAverageVelocity,
    METH_VARARGS, 
    "add doc"},
  { "accumulateExteriorElementPressureIntegrals", 
    cfemIntegralsAccumulateExteriorElementPressureIntegrals,
    METH_VARARGS, 
    "compute average pressure over different exterior element boundaries requires pressure and boundary measure to be zeroed outside call"},
  { "calculateConservationResidual", 
    cfemIntegralsCalculateConservationResidual,
    METH_VARARGS, 
    "calculate discrete mass balance error for a given velocity" },
  { "calculateConservationResidualDG", 
    cfemIntegralsCalculateConservationResidualDG,
    METH_VARARGS, 
    "calculate discrete mass balance error for DG method" },
  {"loadBoundaryFluxIntoGlobalElementBoundaryVelocity",
   cfemIntegralsLoadBoundaryFluxIntoGlobalElementBoundaryVelocity,
   METH_VARARGS,
   "copy flux into boundary velocities with v = updateCoef*v + flux"},
  {"copyGlobalElementBoundaryVelocityToElementBoundary",
   cfemIntegralsCopyGlobalElementBoundaryVelocityToElementBoundary,
   METH_VARARGS,
   "copy global velocity onto element boundary velocities"},
  {"estimate_mt",
   cfemIntegralsEstimate_mt,
   METH_VARARGS,
   "estimate the value of mt at the intial conditions"},
  {"estimate_mt_lowmem",
   cfemIntegralsEstimate_mt_lowmem,
   METH_VARARGS,
   "estimate the value of mt at the intial conditions"},
  {"copyExteriorElementBoundaryValuesFromElementBoundaryValues",
   cfemIntegralsCopyExteriorElementBoundaryValuesFromElementBoundaryValues,
   METH_VARARGS,
   "copy quantity in element boundary quadrature array to one that lives only on exterior element boundaries"},
  {"copyExteriorElementBoundaryValuesToElementBoundaryValues",
   cfemIntegralsCopyExteriorElementBoundaryValuesToElementBoundaryValues,
   METH_VARARGS,
   "copy quantity in element boundary quadrature array from one that lives only on exterior element boundaries"},
  {"copyExteriorElementBoundaryValuesToGlobalElementBoundaryValues",
   cfemIntegralsCopyExteriorElementBoundaryValuesToGlobalElementBoundaryValues,
   METH_VARARGS,
   "copy quantity in global element boundary quadrature array from one that lives only on exterior element boundaries"},
  {"copyExteriorElementBoundaryValuesFromGlobalElementBoundaryValues",
   cfemIntegralsCopyExteriorElementBoundaryValuesFromGlobalElementBoundaryValues,
   METH_VARARGS,
   "copy quantity in global element boundary quadrature array to one that lives only on exterior element boundaries"},
  {"scalarDomainIntegral",
   cfemIntegralsScalarDomainIntegral,
   METH_VARARGS,
   "calculate the integral of a scalar over the domain"},
  {"scalarHeavisideDomainIntegral",
   cfemIntegralsScalarHeavisideDomainIntegral,
   METH_VARARGS,
   "calculate the integral of a scalar over the domain"},
  {"scalarSmoothedHeavisideDomainIntegral",
   cfemIntegralsScalarSmoothedHeavisideDomainIntegral,
   METH_VARARGS,
   "calculate the integral of a scalar over the domain"},
  {"fluxDomainBoundaryIntegral",
   cfemIntegralsFluxDomainBoundaryIntegral,
   METH_VARARGS,
   "calculate the integral of a flux over the boundary of the domain"},
  {"fluxDomainBoundaryIntegralFromVector",
   cfemIntegralsFluxDomainBoundaryIntegralFromVector,
   METH_VARARGS,
   "calculate the integral of a flux over the boundary of the domain"},
 {"computeC0P1InterpolantDGP0",
  cfemIntegralsComputeC0P1InterpolantDGP0,
  METH_VARARGS,
  "compute nodal C0P1 interpolant for DGP0 function "},
 {"computeC0P1InterpolantDGP12",
  cfemIntegralsComputeC0P1InterpolantDGP12,
  METH_VARARGS,
  "compute nodal C0P1 interpolant for DGP1,2 function (Lagrange shape functions)"},
 {"computeC0P1InterpolantNCP1",
  cfemIntegralsComputeC0P1InterpolantNCP1,
  METH_VARARGS,
  "compute nodal C0P1 interpolant for Nonconforming P1 (Crouzeix-Raviart shape functions)"},
 {"checkElementBoundaryAndExteriorElementBoundaryArraysSame",
  cfemIntegralsCheckElementBoundaryAndExteriorElementBoundaryArraysSame,
  METH_VARARGS,
  "check for agreement in  element boundary quadrature array with one that lives only on exterior element boundaries"},
 {"checkGlobalElementBoundaryAndExteriorElementBoundaryArraysSame",
  cfemIntegralsCheckGlobalElementBoundaryAndExteriorElementBoundaryArraysSame,
  METH_VARARGS,
  "check for agreement in  element boundary quadrature array with one that lives only on exterior element boundaries"},
 {"calculateExteriorElementBoundaryStress3D",
  cfemIntegralsCalculateExteriorElementBoundaryStress3D,
  METH_VARARGS,
  "calculate momentum flux  on boundaries for stokes/navier-stokes"},
 {"calculateExteriorElementBoundaryStress2D",
  cfemIntegralsCalculateExteriorElementBoundaryStress2D,
  METH_VARARGS,
  "calculate momentum flux on boundaries for stokes/navier-stokes"},
 {"copyBetweenFreeUnknownsAndGlobalUnknowns",
  cfemIntegralsCopyBetweenFreeUnknownsAndGlobalUnknowns,
  METH_VARARGS,
  "copy over free degrees of freedom to solution"},
  {"updateInteriorElementBoundaryDiffusionAdjoint",
   cfemIntegralsUpdateInteriorElementBoundaryDiffusionAdjoint,
   METH_VARARGS,
   "update the boundary adjoint term for diffusion"},
  {"updateExteriorElementBoundaryDiffusionAdjoint",
   cfemIntegralsUpdateExteriorElementBoundaryDiffusionAdjoint,
   METH_VARARGS,
   "update the boundary adjoint term for diffusion"},
  {"updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense",
   cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense,
   METH_VARARGS,
   "update the Jacobian (dense) from the boundary adjoint term for diffusion (interior boundaries)"},
 {"updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense",
  cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense,
  METH_VARARGS,
  "update the Jacobian (dense) fromt he boundary adjoint term for diffusion (exterior boundaries)"},
  {"updateInteriorElementBoundaryDiffusionAdjoint_sd",
   cfemIntegralsUpdateInteriorElementBoundaryDiffusionAdjoint_sd,
   METH_VARARGS,
   "update the boundary adjoint term for diffusion"},
  {"updateExteriorElementBoundaryDiffusionAdjoint_sd",
   cfemIntegralsUpdateExteriorElementBoundaryDiffusionAdjoint_sd,
   METH_VARARGS,
   "update the boundary adjoint term for diffusion"},
  {"updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense_sd",
   cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_dense_sd,
   METH_VARARGS,
   "update the Jacobian (dense) from the boundary adjoint term for diffusion (interior boundaries)"},
 {"updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense_sd",
  cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_dense_sd,
  METH_VARARGS,
  "update the Jacobian (dense) fromt he boundary adjoint term for diffusion (exterior boundaries)"},
  {"updateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_CSR_sd",
   cfemIntegralsUpdateGlobalJacobianFromInteriorElementBoundaryDiffusionAdjoint_CSR_sd,
   METH_VARARGS,
   "update the Jacobian (dense) from the boundary adjoint term for diffusion (interior boundaries)"},
 {"updateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_CSR_sd",
  cfemIntegralsUpdateGlobalJacobianFromExteriorElementBoundaryDiffusionAdjoint_CSR_sd,
  METH_VARARGS,
  "update the Jacobian (dense) fromt he boundary adjoint term for diffusion (exterior boundaries)"},
 {"update_f_movingDomain",
  cfemIntegralsUpdate_f_movingDomain,
  METH_VARARGS,
  "update the advection term for a moving domain"},
 {"update_f_movingDomain_constantMass",
  cfemIntegralsUpdate_f_movingDomain_constantMass,
  METH_VARARGS,
  "update the advection term for a moving domain when the mass term is 1"},
 {"update_df_movingDomain",
  cfemIntegralsUpdate_df_movingDomain,
  METH_VARARGS,
  "update the derivative with respect to u of the advection term for a moving domain"},
 {"updateStress_weak",
  cfemIntegralsUpdateStress_weak,
  METH_VARARGS,
  "update the weak stress term"},
 {"updateStressJacobian_weak",
  cfemIntegralsUpdateStressJacobian_weak,
  METH_VARARGS,
  "update the Jacobian of the weak stress term"},
 { "projectFromNodalInterpolationConditions",
    cfemIntegralsProjectFromNodalInterpolationConditions,
    METH_VARARGS, 
    "load nodal degrees of freedom from an element based set of interpolation values"},
  { NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initcfemIntegrals(void)
{
  PyObject *m,*d;
  m = Py_InitModule("cfemIntegrals", cfemIntegralsMethods);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
