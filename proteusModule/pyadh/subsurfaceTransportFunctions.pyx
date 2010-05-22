import numpy
cimport numpy
cdef extern from "math.h":
   double fabs(double x)

cdef inline double double_max(double a, double b): return a if a >= b else b
cdef inline double double_min(double a, double b): return a if a <= b else b

ctypedef numpy.double_t DTYPE_t
#numpy.intc_t not in cython's numpy.pxd 
ctypedef int ITYPE_t

def setExteriorElementBoundaryTypes(int nExteriorElementBoundaries_global, 
                                    numpy.ndarray[ITYPE_t,ndim=1] exteriorElementBoundariesArray,
                                    numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryElementsArray,
                                    numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes,
                                    numpy.ndarray[ITYPE_t,ndim=1] exteriorElementBoundaryMaterialTypes):
    cdef int ebNE,ebN,eN
    for ebNE in range(nExteriorElementBoundaries_global):
        ebN = exteriorElementBoundariesArray[ebNE]
        eN  = elementBoundaryElementsArray[ebN,0]
        exteriorElementBoundaryMaterialTypes[ebNE] = elementMaterialTypes[eN]

def setElementBoundariesArray(int nElementBoundaries_global,
                              numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryElementsArray,
                              numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes,
                              numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryMaterialTypes):
    cdef int ebN,eN_left,eN_right
    for ebN in range(nElementBoundaries_global):
        eN_left = elementBoundaryElementsArray[ebN,0]
        eN_right= elementBoundaryElementsArray[ebN,1]
        elementBoundaryMaterialTypes[ebN,0] = elementMaterialTypes[eN_left]
        elementBoundaryMaterialTypes[ebN,1] = elementMaterialTypes[eN_left]
        if eN_right >= 0:
            elementBoundaryMaterialTypes[ebN,1] = elementMaterialTypes[eN_right]

###
def setScalarMaterialFunctionOverElements(numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes,
                                          numpy.ndarray[DTYPE_t,ndim=2] q_vals,
                                          dict material_functions):
    """
    loop over quadrature array and set is material j
    likely little improvement right now without correct typing of material_functions
    """
    cdef int eN,k,material
    for eN in range(q_vals.shape[0]):
        material = elementMaterialTypes[eN]
        for k in range(q_vals.shape[1]):
            q_vals[eN,k] = material_functions[material]

def setVectorMaterialFunctionOverElements(numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes,
                                          numpy.ndarray[DTYPE_t,ndim=3] q_vals,
                                          dict material_functions):
    """
    loop over quadrature array and set \vec f_j assuming element is material j
    """
    cdef int eN,k,material
    for eN in range(q_vals.shape[0]):
        material = elementMaterialTypes[eN]
        for k in range(q_vals.shape[1]):
            q_vals[eN,k,:] = material_functions[material].flat
                                    

def setScalarMaterialFunctionOverElementBoundaries_arithmeticAverage(numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                     numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                     numpy.ndarray[DTYPE_t,ndim=3] ebq_vals,
                                                                     dict material_functions):
    """
    loop over quadrature array and set f = 0.5(f^L_j+f^R_k) assuming element on left 
    is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int eN,ebN,ebN_local,k,material_left,material_right
    
    for eN in range(ebq_vals.shape[0]):
        for ebN_local in range(ebq_vals.shape[1]):
            ebN = elementBoundariesArray[eN,ebN_local]
            material_left = elementBoundaryTypes[ebN,0]
            material_right= elementBoundaryTypes[ebN,1]
            for k in range(ebq_vals.shape[2]):
                ebq_vals[eN,ebN_local,k] = 0.5*(material_functions[material_left]+
                                                material_functions[material_right])

def setSparseTensorMaterialFunctionOverElementBoundaries_harmonicAverage(int nd,
                                                                         numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                         numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                         numpy.ndarray[DTYPE_t,ndim=4] ebq_vals,
                                                                         dict material_functions):
    """
    loop over quadrature array and evaluate function \ten f_{mn} = f^L_{j,mn} f^R_{k,mn}/(f^L_{j,mn}+f^R_{k,mn})
     assuming element on left is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int eN,ebN,ebN_local,k,material_left,material_right,I,J
    cdef double numer,denom

    for eN in range(ebq_vals.shape[0]):
        for ebN_local in range(ebq_vals.shape[1]):
            ebN = elementBoundariesArray[eN,ebN_local]
            material_left = elementBoundaryTypes[ebN,0]
            material_right= elementBoundaryTypes[ebN,1]
            for k in range(ebq_vals.shape[2]):
                for I in range(nd):
                    for J in range(nd):
                        numer = 2.0*material_functions[material_left][I,J]*material_functions[material_right][I,J]
                        denom = material_functions[material_left][I,J] + material_functions[material_right][I,J] + 1.0e-20
                        ebq_vals[eN,ebN_local,k,I*nd+J] = numer/denom

def setScalarMaterialFunctionOverGlobalElementBoundaries_arithmeticAverage(numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                           numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                           numpy.ndarray[DTYPE_t,ndim=2] ebq_global_vals,
                                                                           dict material_functions):
    """
    loop over quadrature array and evaluate function f = 0.5(f^L_j+f^R_k) assuming element on left 
    is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int ebN,material_left,material_right
    
    for ebN in range(ebq_global_vals.shape[0]):
        material_left = elementBoundaryTypes[ebN,0]
        material_right= elementBoundaryTypes[ebN,1]
        for k in range(ebq_global_vals.shape[1]):
            ebq_global_vals[ebN,k] = 0.5*(material_functions[material_left]+
                                          material_functions[material_right])


def setSparseTensorMaterialFunctionOverGlobalElementBoundaries_harmonicAverage(int nd,
                                                                               numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                               numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                               numpy.ndarray[DTYPE_t,ndim=3] ebq_global_vals,
                                                                               dict material_functions):
    """
    loop over quadrature array and evaluate function \ten f_{mn} = f^L_{j,mn}f^R_{k,mn}/(f^L_{j,mn}+f^R_{k,mn})
     assuming element on left is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int ebN,k,material_left,material_right,I,J
    cdef double numer,denom

    for ebN in range(ebq_global_vals.shape[0]):
        material_left = elementBoundaryTypes[ebN,0]
        material_right= elementBoundaryTypes[ebN,1]
        for k in range(ebq_global_vals.shape[1]):
            for I in range(nd):
                for J in range(nd):
                    numer = 2.0*material_functions[material_left][I,J]*material_functions[material_right][I,J]
                    denom = material_functions[material_left][I,J] + material_functions[material_right][I,J] + 1.0e-20
                    ebq_global_vals[ebN,k,I*nd+J] = numer/denom

###
def evaluateScalarMaterialFunctionOverElements(double t,
                                               numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes,
                                               numpy.ndarray[DTYPE_t,ndim=3] x,
                                               numpy.ndarray[DTYPE_t,ndim=2] q_vals,
                                               dict material_functions):
    """
    loop over quadrature array and evaluate function f_j(x,t) assuming element is material j
    likely little improvement right now without correct typing of material_functions
    """
    cdef int eN,k,material
    for eN in range(x.shape[0]):
        material = elementMaterialTypes[eN]
        for k in range(x.shape[1]):
            q_vals[eN,k] = material_functions[material](x[eN,k],t)
                                    
def evaluateVectorMaterialFunctionOverElements(double t,
                                               numpy.ndarray[ITYPE_t,ndim=1] elementMaterialTypes,
                                               numpy.ndarray[DTYPE_t,ndim=3] x,
                                               numpy.ndarray[DTYPE_t,ndim=3] q_vals,
                                               dict material_functions):
    """
    loop over quadrature array and evaluate function \vec f_j(x,t) assuming element is material j
    """
    cdef int eN,k,material
    for eN in range(x.shape[0]):
        material = elementMaterialTypes[eN]
        for k in range(x.shape[1]):
            q_vals[eN,k,:] = material_functions[material](x[eN,k],t).flat
                                    

def evaluateScalarMaterialFunctionOverElementBoundaries_arithmeticAverage(double t,
                                                                          numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                          numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                          numpy.ndarray[DTYPE_t,ndim=4] x,
                                                                          numpy.ndarray[DTYPE_t,ndim=3] ebq_vals,
                                                                          dict material_functions):
    """
    loop over quadrature array and evaluate function f(x,t) = 0.5(f^L_j(x,t)+f^R_k(x,t)) assuming element on left 
    is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int eN,ebN,ebN_local,k,material_left,material_right
    
    for eN in range(x.shape[0]):
        for ebN_local in range(x.shape[1]):
            ebN = elementBoundariesArray[eN,ebN_local]
            material_left = elementBoundaryTypes[ebN,0]
            material_right= elementBoundaryTypes[ebN,1]
            for k in range(x.shape[2]):
                ebq_vals[eN,ebN_local,k] = 0.5*(material_functions[material_left](x[eN,ebN_local,k],t)+
                                                material_functions[material_right](x[eN,ebN_local,k],t))

def evaluateSparseTensorMaterialFunctionOverElementBoundaries_harmonicAverage(int nd,
                                                                              double t,
                                                                              numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                              numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                              numpy.ndarray[DTYPE_t,ndim=4] x,
                                                                              numpy.ndarray[DTYPE_t,ndim=4] ebq_vals,
                                                                              dict material_functions):
    """
    loop over quadrature array and evaluate function \ten f_{mn}(x,t) = f^L_{j,mn}(x,t)f^R_{k,mn}(x,t)/(f^L_{j,mn}(x,t)+f^R_{k,mn})
     assuming element on left is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int eN,ebN,ebN_local,k,material_left,material_right,I,J
    cdef double numer,denom

    for eN in range(x.shape[0]):
        for ebN_local in range(x.shape[1]):
            ebN = elementBoundariesArray[eN,ebN_local]
            material_left = elementBoundaryTypes[ebN,0]
            material_right= elementBoundaryTypes[ebN,1]
            for k in range(x.shape[2]):
                for I in range(nd):
                    for J in range(nd):
                        numer = 2.0*material_functions[material_left](x[eN,ebN_local,k],t)[I,J]*material_functions[material_right](x[eN,ebN_local,k],t)[I,J]
                        denom = material_functions[material_left](x[eN,ebN_local,k],t)[I,J] + material_functions[material_right](x[eN,ebN_local,k],t)[I,J] + 1.0e-20
                        ebq_vals[eN,ebN_local,k,I*nd+J] = numer/denom

def evaluateScalarMaterialFunctionOverGlobalElementBoundaries_arithmeticAverage(double t,
                                                                                numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                                numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                                numpy.ndarray[DTYPE_t,ndim=3] x,
                                                                                numpy.ndarray[DTYPE_t,ndim=2] ebq_global_vals,
                                                                                dict material_functions):
    """
    loop over quadrature array and evaluate function f(x,t) = 0.5(f^L_j(x,t)+f^R_k(x,t)) assuming element on left 
    is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int ebN,material_left,material_right
    
    for ebN in range(x.shape[0]):
        material_left = elementBoundaryTypes[ebN,0]
        material_right= elementBoundaryTypes[ebN,1]
        for k in range(x.shape[1]):
            ebq_global_vals[ebN,k] = 0.5*(material_functions[material_left](x[ebN,k],t)+
                                          material_functions[material_right](x[ebN,k],t))


def evaluateSparseTensorMaterialFunctionOverGlobalElementBoundaries_harmonicAverage(int nd,
                                                                                    double t,
                                                                                    numpy.ndarray[ITYPE_t,ndim=2] elementBoundariesArray,
                                                                                    numpy.ndarray[ITYPE_t,ndim=2] elementBoundaryTypes,
                                                                                    numpy.ndarray[DTYPE_t,ndim=4] x,
                                                                                    numpy.ndarray[DTYPE_t,ndim=3] ebq_global_vals,
                                                                                    dict material_functions):
    """
    loop over quadrature array and evaluate function \ten f_{mn}(x,t) = f^L_{j,mn}(x,t)f^R_{k,mn}(x,t)/(f^L_{j,mn}(x,t)+f^R_{k,mn})
     assuming element on left is material j and element on right is material k

    likely little improvement right now without correct typing of material_functions
 
    """
    cdef int ebN,k,material_left,material_right,I,J
    cdef double numer,denom

    for ebN in range(x.shape[0]):
        material_left = elementBoundaryTypes[ebN,0]
        material_right= elementBoundaryTypes[ebN,1]
        for k in range(x.shape[1]):
            for I in range(nd):
                for J in range(nd):
                    numer = 2.0*material_functions[material_left](x[ebN,k],t)[I,J]*material_functions[material_right](x[ebN,k],t)[I,J]
                    denom = material_functions[material_left](x[ebN,k],t)[I,J] + material_functions[material_right](x[ebN,k],t)[I,J] + 1.0e-20
                    ebq_global_vals[ebN,k,I*nd+J] = numer/denom

