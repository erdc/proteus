#include "Python.h"
#include "numpy/arrayobject.h"

/** \file triangleWrappersModule.c
    \defgroup triangleWrappers triangleWrappers
    \brief Python interface to triangle
    @{
*/

/***********************************************************************
  try to write an interface to Triangle in python. I will follow the
  ellipt2d example as much as possible, but I will setup a Numeric
  interface for the python data types coming in

  
**********************************************************************/
#define REAL double
#define _NDIM 2
/*Triangle package's interface: triangulateio */
#include PROTEUS_TRIANGLE_H



/*taken directly from ellipt2d*/
typedef void (*destr)(void *);
void destroy_triangulateio(struct triangulateio *object){
#if defined(Py_DEBUG) || defined(DEBUG)
  printf("now destroying triangulateio\n");
#endif
  /*mwf now set these to NULL after freeing*/
  if( object->pointlist             ) 
    {
      free( object->pointlist             );
      object->pointlist = NULL;
    }
  if( object->pointattributelist    ) 
    {
      free( object->pointattributelist    );
      object->pointattributelist = NULL;
    }
  if( object->pointmarkerlist       ) 
    {
      free( object->pointmarkerlist       );
      object->pointmarkerlist = NULL;
    }
  if( object->trianglelist          ) 
    {
      free( object->trianglelist          );
      object->trianglelist = NULL;
    }
  if( object->triangleattributelist ) 
    {
      free( object->triangleattributelist );
      object->triangleattributelist = NULL;
    }
  if( object->trianglearealist      ) 
    {
      free( object->trianglearealist      );
      object->trianglearealist = NULL;
    }
  if( object->neighborlist          ) 
    {
      free( object->neighborlist          );
      object->neighborlist = NULL;
    }
  if( object->segmentlist           ) 
    {
      free( object->segmentlist           );
      object->segmentlist = NULL;
    }
  if( object->segmentmarkerlist     ) 
    {
      free( object->segmentmarkerlist     );
      object->segmentmarkerlist = NULL;
    }
  /*mwf holelist and regionlist are shallow copies after
        triangulate is called, so I need to make sure they
        get freed only once somehow*/
  if( object->holelist && object->holelistmemflag == 1) 
    {
      free( object->holelist              );
      object->holelist = NULL;
    }
  if( object->regionlist && object->regionlistmemflag == 1) 
    {
      free( object->regionlist            );
      object->regionlist = NULL;
    }
  if( object->edgelist              ) 
    {
      free( object->edgelist              );
      object->edgelist = NULL;
    }
  if( object->edgemarkerlist        ) 
    {
      free( object->edgemarkerlist        );
      object->edgemarkerlist = NULL;
    }
  if( object->normlist              ) 
    {
      free( object->normlist              );
      object->normlist = NULL;
    }
  free(object);
  object = NULL;
}


/***********************************************************************
  Automate the conversion of certain array types since this is 
  essentially what most of the operations are

  triArray should point to the array in triangulateio object
    that we are setting
  triSize[0,1] are the corresponding dimensions in the
    triangulateio object
  triSizeSet is the corresponding dimension that should be set in 
    triangulateio object
    -1 means set both, 0, 1 is the dimension to set, anthing else means ignore
  pyArray is the Python object that contains the data to be copied
  forceResize = 1 means always redimension the triangulateio data member
 **********************************************************************/
static PyObject *
copyToTriangulateioDoubleArray1or2dim(double ** ptriArray, 
				      int *triSize0, int *triSize1,
				      PyObject *pyArray,
				      int forceResize, int triSizeSet)
{
  PyArrayObject *array;
  double *araw;

  int i,j;
  int dim0,dim1;
  int minDim,maxDim;

  assert(-2 <= triSizeSet && triSizeSet <= 1);
  minDim = 1;
  maxDim = 2;
  array = (PyArrayObject *) PyArray_ContiguousFromObject(pyArray,
							 PyArray_DOUBLE,
							 minDim,maxDim);

 if (!array)
    {
      PyErr_Format(PyExc_ValueError,"copyToTriDoubleArray incorrect dim:%d:%d \n",
		   minDim,maxDim);
      return NULL;
    }
 dim0 = array->dimensions[0];
 dim1 = 1;
 if (array->nd > 1)
   dim1 = array->dimensions[1];
 if (dim0 > 0)
   {
     if (*triSize0 != dim0 || *triSize1 != dim1 || forceResize > 0)
       {
	 if (*ptriArray) free(*ptriArray);
	 *ptriArray = (REAL*)malloc(dim0*dim1*sizeof(REAL));
       }
     if (triSizeSet == 1)
       *triSize1 = dim1;
     else if (triSizeSet == 0)
       *triSize0 = dim0;
     else if (triSizeSet == -1)
       {
	 *triSize1 = dim1; *triSize0 = dim0;
       }
   }

 /*now copy over data*/
 araw = (double *) array->data;
#ifdef DEBUG_UTIL
 printf("copyToTriDoubleArray dim0= %d dim1= %d triSizeSet= %d ptriArray= %x\n",
	dim0,dim1,triSizeSet,*ptriArray);
#endif

 for (i = 0; i < dim0; i++)
   {
     for (j = 0; j < dim1; j++)
       {
	 (*ptriArray)[i*dim1 + j] = araw[i*dim1 + j];
#ifdef DEBUG_UTIL
	 printf("ptriArray[%d,%d]= %g \n",i,j,(*ptriArray)[i*dim1+j]);
#endif
       }
   }

  /*free up arrays created*/
  Py_DECREF(array);
  return Py_BuildValue("");

}
/***********************************************************************
  Automate the conversion of certain array types since this is 
  essentially what most of the operations are
  This version is for integer arrays

  ptriArray should point to the array in triangulateio object
    that we are setting
  triSize[0,1] are the corresponding dimensions in the
    triangulateio object
  triSizeSet is the corresponding dimension that should be set in 
    triangulateio object
    -1 means set both, 0, 1 is the dimension to set, anthing else means ignore

  pyArray is the Python object that contains the data to be copied
  forceResize = 1 means always redimension the triangulateio data member
 **********************************************************************/
static PyObject *
copyToTriangulateioIntegerArray1or2dim(int ** ptriArray, 
				       int *triSize0, int *triSize1,
				       PyObject *pyArray,
				       int forceResize, int triSizeSet)
{
  PyArrayObject *array;
  int *araw;

  int i,j;
  int dim0,dim1;
  int minDim,maxDim;

  assert(-2 <= triSizeSet && triSizeSet <= 1);
  minDim = 1;
  maxDim = 2;
  array = (PyArrayObject *) PyArray_ContiguousFromObject(pyArray,
							 PyArray_INT,
							 minDim,maxDim);

 if (!array)
    {
      PyErr_Format(PyExc_ValueError,"copyToTriIntArray incorrect dim:%d:%d \n",
		   minDim,maxDim);
      return NULL;
    }
 dim0 = array->dimensions[0];
 dim1 = 1;
 if (array->nd > 1)
   dim1 = array->dimensions[1];
 if (dim0 > 0)
   {
     if (*triSize0 != dim0 || *triSize1 != dim1 || forceResize > 0)
       {
	 if (*ptriArray) free(*ptriArray);
	 *ptriArray = (int*)malloc(dim0*dim1*sizeof(int));
       }
     if (triSizeSet == 1)
       *triSize1 = dim1;
     else if (triSizeSet == 0)
       *triSize0 = dim0;
     else if (triSizeSet == -1)
       {
	 *triSize1 = dim1; *triSize0 = dim0;
       }
   }

 /*now copy over data*/
 araw = (int *) array->data;
#ifdef DEBUG_UTIL
	dim0,dim1,triSizeSet,*ptriArray);
#endif

 for (i = 0; i < dim0; i++)
   {
     for (j = 0; j < dim1; j++)
       {
	 (*ptriArray)[i*dim1 + j] = araw[i*dim1 + j];
#ifdef DEBUG_UTIL
	 printf("ptriArray[%d,%d]= %d \n",i,j,(*ptriArray)[i*dim1+j]);
#endif
       }
   }

  /*free up arrays created*/
  Py_DECREF(array);
  return Py_BuildValue("");

}

/***********************************************************************
  Automate getting copies of arrays from triangulateio to send to python
  triArray should point to the array in triangulateio object
    that we are setting
  triSize[0,1] are the corresponding dimensions in the
    triangulateio object
  deepCopy = 1 means make a real copy
 **********************************************************************/
static PyObject * 
getCopyOfTriangulateioDoubleArray1or2dim(double *triArray,
					 int triSize0, int triSize1,
					 int deepCopy)
{
  PyArrayObject *array;
  double *araw;

  int i,j,nd,dims[2];

  nd = 0;
  if (triSize0 > 0)
    nd = 1;
  if (triSize1 > 0)
    nd = 2;
  
  if (nd > 0 && triArray)
    {
      assert(triArray);

      dims[0] = triSize0;
      dims[1] = 1;
      if (nd > 1)
	dims[1] = triSize1;

      if (deepCopy == 1)
	{
	  array = (PyArrayObject *) PyArray_FromDims(nd,dims,PyArray_DOUBLE);
	  if (array == NULL)
	    {
	      printf("getCopyOfDoubleArray1or2dim failed to create ou array dims=%dx%d\n",
		     dims[0],dims[1]);
	      return NULL;
	    }
	  araw  = (double *) array->data;
	  for (i = 0; i < dims[0]; i++)
	    {
	      for (j = 0; j < dims[1]; j++)
		{
		  araw[i*dims[1]+j] = triArray[i*dims[1]+j];
		}
	    }
	}/*end deep copy*/
      else
	{
	  array = (PyArrayObject *) PyArray_FromDimsAndData(nd,dims,PyArray_DOUBLE,
							    (char *) triArray);

	  if (array == NULL)
	    {
	      printf("getCopyOfDoubleArray1or2dim failed to create ou array dims=%dx%d\n",
		     dims[0],dims[1]);
	      return NULL;
	    }
	}/*end shallow copy*/
#ifdef DEBUG_UTIL
      araw  = (double *) array->data;
      /*mwf debug*/
      printf("in getCopyOfDoubleArray1or2dim triArray= %x dim=%dx%d \n",triArray,
	     dims[0],dims[1]);
      for (i=0; i < dims[0]; i++)
	{
	  for (j=0; j < dims[1]; j++)
	    {
	      printf("array[%d,%d]= %g\n",i,j,araw[i*dims[1]+j]);

	    }
	}
      /*mwf debug end*/
#endif
      return PyArray_Return(array);
    }/* end nd > 0*/
  return Py_BuildValue(""); 

}
/***********************************************************************
  Automate getting copies of arrays from triangulateio to send to python
  triArray should point to the array in triangulateio object
    that we are setting
  triSize[0,1] are the corresponding dimensions in the
    triangulateio object
  deepCopy = 1 means make a real copy
 **********************************************************************/
static PyObject * 
getCopyOfTriangulateioIntegerArray1or2dim(int *triArray,
					  int triSize0, int triSize1,
					  int deepCopy)
{
  PyArrayObject *array;
  int *araw;

  int i,j,nd,dims[2];

  nd = 0;
  if (triSize0 > 0)
    nd = 1;
  if (triSize1 > 0)
    nd = 2;
  
  if (nd > 0 && triArray)
    {
      assert(triArray);

      dims[0] = triSize0;
      dims[1] = 1;
      if (nd > 1)
	dims[1] = triSize1;

      if (deepCopy == 1)
	{
	  array = (PyArrayObject *) PyArray_FromDims(nd,dims,PyArray_INT);
	  if (array == NULL)
	    {
	      printf("getCopyOfIntArray1or2dim failed to create ou array dims=%dx%d\n",
		     dims[0],dims[1]);
	      return NULL;
	    }
	  araw  = (int *) array->data;
	  for (i = 0; i < dims[0]; i++)
	    {
	      for (j = 0; j < dims[1]; j++)
		{
		  araw[i*dims[1]+j] = triArray[i*dims[1]+j];
		}
	    }
	}/*end deep copy*/
      else
	{
	  array = (PyArrayObject *) PyArray_FromDimsAndData(nd,dims,PyArray_INT,
							    (char *) triArray);

	  if (array == NULL)
	    {
	      printf("getCopyOfIntArray1or2dim failed to create ou array dims=%dx%d\n",
		     dims[0],dims[1]);
	      return NULL;
	    }
	}/*end shallow copy*/
#ifdef DEBUG_UTIL
      araw  = (int *) array->data;
      /*mwf debug*/
      printf("in getCopyOfDoubleArray1or2dim triArray= %x dim=%dx%d \n",triArray,
	     dims[0],dims[1]);
      for (i=0; i < dims[0]; i++)
	{
	  for (j=0; j < dims[1]; j++)
	    {
	      printf("array[%d,%d]= %d\n",i,j,araw[i*dims[1]+j]);

	    }
	}
      /*mwf debug end*/
#endif
      return PyArray_Return(array);
    }/* end nd > 0*/
  return Py_BuildValue(""); 

}

/***************************************************
 just write some things out
 **************************************************/
void simplePrintReport(struct triangulateio* io, int markers, int reporttriangles,
		       int reportneighbors, int reportsegments,
		       int reportedges, int reportnorms)
{
  int i, j;
  assert(io);

  for (i = 0; i < io->numberofpoints; i++) {
    printf("Point %4d:", i);
    for (j = 0; j < 2; j++) {
      printf("  %.6g", io->pointlist[i * 2 + j]);
    }
    if (io->numberofpointattributes > 0) {
      printf("   attributes");
    }
    for (j = 0; j < io->numberofpointattributes; j++) {
      printf("  %.6g",
             io->pointattributelist[i * io->numberofpointattributes + j]);
    }
    if (markers && io->pointmarkerlist) {
      printf("   marker %d\n", io->pointmarkerlist[i]);
    } else {
      printf("\n");
    }
  }
  printf("\n");

  if (reporttriangles || reportneighbors) {
    for (i = 0; i < io->numberoftriangles; i++) {
      if (reporttriangles && io->trianglelist) {
        printf("Triangle %4d points:", i);
        for (j = 0; j < io->numberofcorners; j++) {
          printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
        }
        if (io->numberoftriangleattributes > 0) {
          printf("   attributes");
        }
        for (j = 0; j < io->numberoftriangleattributes; j++) {
          printf("  %.6g", io->triangleattributelist[i *
                                         io->numberoftriangleattributes + j]);
        }
        printf("\n");
      }
      if (reportneighbors && io->neighborlist) {
        printf("Triangle %4d neighbors:", i);
        for (j = 0; j < 3; j++) {
          printf("  %4d", io->neighborlist[i * 3 + j]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportsegments && io->segmentlist) {
    for (i = 0; i < io->numberofsegments; i++) {
      printf("Segment %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->segmentlist[i * 2 + j]);
      }
      if (markers && io->segmentmarkerlist) {
        printf("   marker %d\n", io->segmentmarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

  if (reportedges && io->edgelist) {
    for (i = 0; i < io->numberofedges; i++) {
      printf("Edge %4d points:", i);
      for (j = 0; j < 2; j++) {
        printf("  %4d", io->edgelist[i * 2 + j]);
      }
      if (reportnorms && io->normlist && (io->edgelist[i * 2 + 1] == -1)) {
        for (j = 0; j < 2; j++) {
          printf("  %.6g", io->normlist[i * 2 + j]);
        }
      }
      if (markers && io->edgemarkerlist) {
        printf("   marker %d\n", io->edgemarkerlist[i]);
      } else {
        printf("\n");
      }
    }
    printf("\n");
  }

}
/***********************************************************************
  end utility routines
  **********************************************************************/


/***************************************************
  create a new triangulateio object
**************************************************/
static PyObject *
triangulate_NEW(PyObject *self, PyObject *args){
  PyObject *address, *result;
  struct triangulateio *object;

  object = (struct triangulateio*)malloc(sizeof(struct triangulateio));

  object->pointlist             = NULL;
  object->pointattributelist    = NULL; /* In / out */
  object->pointmarkerlist       = NULL; /* In / out */
  object->numberofpoints         = 0;    /* In / out */
  object->numberofpointattributes= 0;    /* In / out */

  object->trianglelist          = NULL; /* In / out */
  object->triangleattributelist = NULL; /* In / out */
  object->trianglearealist      = NULL; /* In only */
  object->neighborlist          = NULL; /* Out only */
  object->numberoftriangles     = 0;    /* In / out */
  object->numberofcorners       = 0;    /* In / out */
  object->numberoftriangleattributes = 0;    /* In / out */

  object->segmentlist           = NULL; /* In / out */
  object->segmentmarkerlist     = NULL; /* In / out */
  object->numberofsegments      = 0;    /* In / out */

  object->holelist              = NULL; /* In / pointer to array copied out */
  object->numberofholes         = 0;    /* In / copied out */

  object->regionlist            = NULL; /* In / pointer to array copied out */
  object->numberofregions       = 0;    /* In / copied out */

  object->edgelist              = NULL; /* Out only */
  object->edgemarkerlist        = NULL; /* Not used with Voronoi diagram; out only */
  object->normlist              = NULL; /* Used only with Voronoi diagram; out only */
  object->numberofedges         = 0;    /* Out only */

  /*mwf add flags to keep track of holelist and regionlist since they are shallow
    copies in triangle library
    value of 1 means that this object owns memory for its holelist or regionlist
    That would be the case for an input triangulation in the basic examples
   */
  object->holelistmemflag  = 0; /*does not own memory*/
  object->regionlistmemflag= 0; /*does not own memory*/
  /* return opaque handle */

  address = PyCObject_FromVoidPtr(object, (destr) destroy_triangulateio);
  result = Py_BuildValue("O", address);
  Py_DECREF(address);
  return result;
}
/***************************************************
  interface to the triangulate function in triangle library
  
  triangulate(flags,trin,triout,vorout)

  takes an existing triangulateio object in and processes based on switches
    specified to produce an output triangulateio object and a Voronoi diagram
    output. The input and output objects have to be allocated already
 
**************************************************/
static PyObject *
triangulate_APPLY_TRIANGULATE(PyObject *self, PyObject *args)
{
  PyObject *address0,*address1,*address2, *flagsin;
  struct triangulateio *in,*out,*vor;
  char *flags;
  if(!PyArg_ParseTuple(args,(char *)"OOOO", 
		       &flagsin,&address0,&address1,&address2)){ 
    return NULL;
  }
  
  if(!PyString_Check(flagsin)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! String required (flags).");
    return NULL;
  } 

  if(!PyCObject_Check(address0)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 2nd argument! CObject required (triangulateio handle).");
    return NULL;
  } 
  if(!PyCObject_Check(address1)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 3rd argument! CObject required (triangulateio handle).");
    return NULL;
  }    
  if(!PyCObject_Check(address2)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 4th argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  flags = PyString_AsString(flagsin);
#ifdef DEBUG
  printf("ApplyTriangulate flags= %s\n",flags);
#endif

  in  = (struct triangulateio*)PyCObject_AsVoidPtr(address0);  
  out = (struct triangulateio*)PyCObject_AsVoidPtr(address1);  
  vor = (struct triangulateio*)PyCObject_AsVoidPtr(address2);

  triangulate(flags,in,out,vor);

#ifdef DEBUG
  printf("ApplyTriangulate done\n");
#endif

  /*hopefully the data have been modified appropriately now*/
  return Py_BuildValue("");
}
/***************************************************
  interface to the triangulate function in triangle library
  
  triangulate(flags,trin,triout,NULL)

  takes an existing triangulateio object in and processes based on switches
    specified to produce an output triangulateio object but not a Voronoi diagram
    output. The input and output objects have to be allocated already
 
**************************************************/
static PyObject *
triangulate_APPLY_TRIANGULATE_NO_VORONOI(PyObject *self, PyObject *args)
{
  PyObject *address0,*address1, *flagsin;
  struct triangulateio *in,*out,*vor;
  char *flags;
  if(!PyArg_ParseTuple(args,(char *)"OOO", 
		       &flagsin,&address0,&address1)){ 
    return NULL;
  }
  
  if(!PyString_Check(flagsin)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! String required (flags).");
    return NULL;
  } 

  if(!PyCObject_Check(address0)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 2nd argument! CObject required (triangulateio handle).");
    return NULL;
  } 
  if(!PyCObject_Check(address1)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 3rd argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  flags = PyString_AsString(flagsin);
#ifdef DEBUG
  printf("ApplyTriangulate flags= %s\n",flags);
#endif

  in  = (struct triangulateio *)PyCObject_AsVoidPtr(address0);  
  out = (struct triangulateio *)PyCObject_AsVoidPtr(address1);  
  vor = (struct triangulateio *) NULL;

  triangulate(flags,in,out,vor);

#ifdef DEBUG
  printf("ApplyTriangulate done\n");
#endif

  /*hopefully the data have been modified appropriately now*/
  return Py_BuildValue("");
}

/***************************************************
   just get basic info about the mesh:
     numberofpoints
     numberofpointattributes
     numberoftriangles
     numberofcorners
     numberoftriangleattributes
     numberofsegments
     numberofholes
     numberofregions
     numberofedges
 **************************************************/
static PyObject *
triangulate_GET_INFO(PyObject *self, PyObject *args)
{
  PyObject *address, *rval;
  struct triangulateio *object;
  if(!PyArg_ParseTuple(args,(char *)"O", 
		       &address)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  

  rval = Py_BuildValue("(iiiiiiiii)",object->numberofpoints,
		       object->numberofpointattributes,
		       object->numberoftriangles,
		       object->numberofcorners,
		       object->numberoftriangleattributes,
		       object->numberofsegments,
		       object->numberofholes,
		       object->numberofregions,
		       object->numberofedges);
  return rval;
}
/***************************************************
 write out some info for the triangulate objects

 start by just writing with printf
 **************************************************/
static PyObject *
triangulate_PRINT_REPORT(PyObject *self, PyObject *args)
{
  PyObject *address;
  struct triangulateio *object;
  int markers,reporttriangles,reportneighbors,reportsegments,
    reportedges,reportnorms;

  if(!PyArg_ParseTuple(args,(char *)"O", 
		       &address)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    
  object = (struct triangulateio *) (struct triangulateio *)PyCObject_AsVoidPtr(address);  

  /*print everything for now*/
  markers=1; reporttriangles=1; reportneighbors=1; reportsegments=1;
  reportedges=1; reportnorms=1;
  /*markers=0; reporttriangles=0; reportneighbors=0; reportsegments=0;
    reportedges=0; reportnorms=0;*/

  simplePrintReport(object,
		    markers,reporttriangles,reportneighbors,
		    reportsegments,reportedges,reportnorms);

  return Py_BuildValue("");
}
/***************************************************
   load points into data array, assume that it is
   a 2 dimensional array : Npoints x SpaceDim
 **************************************************/
static PyObject *
triangulate_SET_POINTS(PyObject *self, PyObject *args){
  PyObject *address, *valsin0, *rval0;
  struct triangulateio *object;
  double **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  int spaceDim;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin0)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->pointlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofpoints;
  spaceDim = 2;
  dim1 = &spaceDim;
  /*which dimension do I need to update */
  setDim = 0; /*number points in domain*/
  forceResize = 0; /*go ahead and require resizing of data arrays*/
  rval0 = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
						valsin0,
						forceResize,setDim);


  return rval0;
}
/***************************************************
   load points and markers into data array, assume points is
   a 2 dimensional array : Npoints x SpaceDim
   markers is Npoints x 1
 **************************************************/
static PyObject *
triangulate_SET_POINTS_AND_MARKERS(PyObject *self, PyObject *args){
  PyObject *address, *valsin0, *valsin1, *rval0, *rval1;
  struct triangulateio *object;
  double **ptriArray;
  int **ptriArrayI;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  int spaceDim,markerDim;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OOO", 
		       &address, &valsin0, &valsin1)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->pointlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofpoints;
  spaceDim = 2;
  dim1 = &spaceDim;
  /*which dimension do I need to update */
  setDim = 0; /*number points in domain*/
  forceResize = 0; /*go ahead and require resizing of data arrays*/
  rval0 = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
						valsin0,
						forceResize,setDim);

  /*now repeat for marker list*/
  /*pick out the array that I need to set*/
  ptriArrayI = &(object->pointmarkerlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofpoints;
  markerDim = 1;
  dim1 = &markerDim;
  /*which dimension do I need to update */
  setDim = -2; /*don'r resize anything*/
  forceResize = 1; /*have to reset data arrays here*/
  rval1 = copyToTriangulateioIntegerArray1or2dim(ptriArrayI,dim0,dim1,
						 valsin1,
						 forceResize,setDim);

  return rval1;
}

/***************************************************
   load pointsmarkers into data array, assume pointmarkers
   is Npoints x 1
 **************************************************/
static PyObject *
triangulate_SET_POINTMARKERS(PyObject *self, PyObject *args){
  PyObject *address, *valsin1, *rval1;
  struct triangulateio *object;
  int **ptriArrayI;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  int markerDim;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin1)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out marker list*/
  /*pick out the array that I need to set*/
  ptriArrayI = &(object->pointmarkerlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofpoints;
  assert(dim0 > 0);
  markerDim = 1;
  dim1 = &markerDim;
  /*which dimension do I need to update */
  setDim = -2; /*don'r resize anything*/
  forceResize = 1; /*have to reset data arrays here*/
  rval1 = copyToTriangulateioIntegerArray1or2dim(ptriArrayI,dim0,dim1,
						 valsin1,
						 forceResize,setDim);

  return rval1;
}
/***************************************************
   get points from triangulateio data structure
   as a 2 dimensional array : Npoints x SpaceDim
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_POINTS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->pointlist;   /*array to copy*/
  nd  = 2;                        /*pointslist is 2d array*/
  dim0= object->numberofpoints;   /*dimensions of array to copy*/
  dim1= _NDIM;                    /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get points from triangulateio data structure
   as a 2 dimensional array : Npoints x SpaceDim
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_POINTS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->pointlist;   /*array to copy*/
  nd  = 2;                        /*pointslist is 2d array*/
  dim0= object->numberofpoints;   /*dimensions of array to copy*/
  dim1= _NDIM;                    /*dimensions of array to copy*/
  deepCopy = 1;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get pointmarkerlist from triangulateio data structure
   as a 2 dimensional array : Nsegments x 1
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_POINTMARKERS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->pointmarkerlist;    /*array to copy*/
  nd  = 1;                                 /*trianglelist is 2d array*/
  dim0= object->numberofpoints;          /*dimensions of array to copy*/
  dim1= 0;                                 /*dimensions of array to copy*/
  deepCopy = 0;                            /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get pointmarkerlist from triangulateio data structure
   as a 2 dimensional array : Npoints x 1
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_POINTMARKERS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->pointmarkerlist;    /*array to copy*/
  nd  = 1;                                 /*trianglelist is 2d array*/
  dim0= object->numberofpoints;          /*dimensions of array to copy*/
  dim1= 0;                                 /*dimensions of array to copy*/
  deepCopy = 0;                            /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   load point attributes into data array, assume that it is
   a 2 dimensional array : Npoints x Nattributes
 **************************************************/
static PyObject *
triangulate_SET_POINT_ATTRIBS(PyObject *self, PyObject *args){
  PyObject *address, *valsin, *rval;
  struct triangulateio *object;
  double **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->pointattributelist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofpoints;
  dim1 = &object->numberofpointattributes;
  /*which dimension do I need to update */
  setDim = 1; /*number of attributes per triangle*/
  forceResize = 1; /*go ahead and require resizing of data arrays*/
  rval = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
					       valsin,
					       forceResize,setDim);

  return rval;
}
/***************************************************
   get point attributes from triangulateio data structure
   as a 2 dimensional array : Npoitns x Nattrib
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_POINT_ATTRIBS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->pointattributelist;            /*array to copy*/
  nd  = 2;                                   /*trianglelist is 2d array*/
  dim0= object->numberofpoints;              /*dimensions of array to copy*/
  dim1= object->numberofpointattributes;     /*dimensions of array to copy*/
  deepCopy = 0;                              /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   get point attributes from triangulateio data structure
   as a 2 dimensional array : Npoitns x Nattrib
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_POINT_ATTRIBS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->pointattributelist;            /*array to copy*/
  nd  = 2;                                   /*trianglelist is 2d array*/
  dim0= object->numberofpoints;              /*dimensions of array to copy*/
  dim1= object->numberofpointattributes;     /*dimensions of array to copy*/
  deepCopy = 1;                              /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}



/***************************************************
   load triangle info into data array, assume that it is
   a 2 dimensional array : Ntri x Ncorner (3 or 6)
 **************************************************/
static PyObject *
triangulate_SET_TRIANGLES(PyObject *self, PyObject *args){
  PyObject *address, *valsin, *rval;
  struct triangulateio *object;
  int **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  /*arguments triangulateioHandle,triangles*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->trianglelist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberoftriangles;
  dim1 = &object->numberofcorners;
  /*which dimension do I need to update */
  setDim = -1; /*set both dimensions*/
  forceResize = 0; /*go ahead and require resizing of data arrays*/
  rval  = copyToTriangulateioIntegerArray1or2dim(ptriArray,dim0,dim1,
						 valsin,
						 forceResize,setDim);

  return rval;
}



/***************************************************
   get triangles from triangulateio data structure
   as a 2 dimensional array : Npoints x Ncorners
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_TRIANGLES(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->trianglelist;   /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberoftriangles;   /*dimensions of array to copy*/
  dim1= object->numberofcorners;     /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get triangles from triangulateio data structure
   as a 2 dimensional array : Npoints x Ncorners
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_TRIANGLES_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->trianglelist;   /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberoftriangles;   /*dimensions of array to copy*/
  dim1= object->numberofcorners;     /*dimensions of array to copy*/
  deepCopy = 1;                      /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}



/***************************************************
   load triangle attributes into data array, assume that it is
   a 2 dimensional array : Nelems x Nattributes
 **************************************************/
static PyObject *
triangulate_SET_TRIANGLE_ATTRIBS(PyObject *self, PyObject *args){
  PyObject *address, *valsin, *rval;
  struct triangulateio *object;
  double **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->triangleattributelist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberoftriangles;
  dim1 = &object->numberoftriangleattributes;
  /*which dimension do I need to update */
  setDim = 1; /*number of attributes per triangle*/
  forceResize = 1; /*go ahead and require resizing of data arrays*/
  rval = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
					       valsin,
					       forceResize,setDim);

  return rval;
}
/***************************************************
   get triangle attributes from triangulateio data structure
   as a 2 dimensional array : Ntriangles x Nattrib
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_TRIANGLE_ATTRIBS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->triangleattributelist;  /*array to copy*/
  nd  = 2;                                   /*trianglelist is 2d array*/
  dim0= object->numberoftriangles;           /*dimensions of array to copy*/
  dim1= object->numberoftriangleattributes;  /*dimensions of array to copy*/
  deepCopy = 0;                              /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get triangle attributes from triangulateio data structure
   as a 2 dimensional array : Ntriangles x Nattrib
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_TRIANGLE_ATTRIBS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->triangleattributelist;  /*array to copy*/
  nd  = 2;                                   /*trianglelist is 2d array*/
  dim0= object->numberoftriangles;           /*dimensions of array to copy*/
  dim1= object->numberoftriangleattributes;   /*dimensions of array to copy*/
  deepCopy = 1;                              /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   load triangle area contraints into data array, assume that it is
   a 1 dimensional array : Nelems x 1
 **************************************************/
static PyObject *
triangulate_SET_TRIANGLE_AREAS(PyObject *self, PyObject *args){
  PyObject *address, *valsin, *rval;
  struct triangulateio *object;
  double **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int scalarDim;
  int forceResize;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->trianglearealist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberoftriangles;
  scalarDim = 1;
  dim1 = &scalarDim;
  /*which dimension do I need to update */
  setDim = -2; /*number of attributes per triangle*/
  forceResize = 1; /*go ahead and require resizing of data arrays*/
  rval = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
					       valsin,
					       forceResize,setDim);

  return rval;
}

/***************************************************
   get triangle neighborlist from triangulateio data structure
   as a 2 dimensional array : Nelems x 3
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_NEIGHBORS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->neighborlist;   /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberoftriangles;   /*dimensions of array to copy*/
  dim1= 3;                           /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get triangle neighborlist from triangulateio data structure
   as a 2 dimensional array : Nelems x 3
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_NEIGHBORS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->neighborlist;   /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberoftriangles;   /*dimensions of array to copy*/
  dim1= 3;                           /*dimensions of array to copy*/
  deepCopy = 1;                      /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   load segment info into data array, assume that it is
   a 2 dimensional array : Nsegments x 2
 **************************************************/
static PyObject *
triangulate_SET_SEGMENTS(PyObject *self, PyObject *args){
  PyObject *address, *valsin0, *rval0;
  struct triangulateio *object;
  int **ptriArray;
  int *dim0,*dim1;
  int setDim; 
  int forceResize;
  int spaceDim;
  /*arguments triangulateioHandle,segments,segmarks*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin0)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->segmentlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofsegments;
  spaceDim = 2;
  dim1 = &spaceDim;
  /*which dimension do I need to update */
  setDim = 0; /*number segments in domain*/
  forceResize = 0; /*go ahead and require resizing of data arrays*/
  rval0 = copyToTriangulateioIntegerArray1or2dim(ptriArray,dim0,dim1,
						 valsin0,
						 forceResize,setDim);

  return rval0;
}
/***************************************************
   load segment info into data array, assume that it is
   a 2 dimensional array : Nsegments x 2
   also load in segment markers : Nsegments x 1
 **************************************************/
static PyObject *
triangulate_SET_SEGMENTS_AND_MARKERS(PyObject *self, PyObject *args){
  PyObject *address, *valsin0, *valsin1, *rval0, *rval1;
  struct triangulateio *object;
  int **ptriArray;
  int **ptriArrayI;
  int *dim0, *dim1;
  int setDim; 
  int forceResize;
  int spaceDim,markerDim;
  /*arguments triangulateioHandle,segments,segmarks*/
  if(!PyArg_ParseTuple(args,(char *)"OOO", 
		       &address, &valsin0, &valsin1)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->segmentlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofsegments;
  spaceDim = 2;
  dim1 = &spaceDim;
  /*which dimension do I need to update */
  setDim = 0; /*number segments in domain*/
  forceResize = 0; /*go ahead and require resizing of data arrays*/
  rval0 = copyToTriangulateioIntegerArray1or2dim(ptriArray,dim0,dim1,
						 valsin0,
						 forceResize,setDim);

  /*now repeat for marker list*/
  /*pick out the array that I need to set*/
  ptriArrayI = &(object->segmentmarkerlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofsegments;
  markerDim = 1;
  dim1 = &markerDim;
  /*which dimension do I need to update */
  setDim = -2; /*don't resize anything*/
  forceResize = 1; /*have to reset data arrays here*/
  rval1 = copyToTriangulateioIntegerArray1or2dim(ptriArrayI,dim0,dim1,
						 valsin1,
						 forceResize,setDim);

  return rval1;
}
/***************************************************
   load in segment markers : Nsegments x 1
 **************************************************/
static PyObject *
triangulate_SET_SEGMENTMARKERS(PyObject *self, PyObject *args){
  PyObject *address, *valsin1, *rval1;
  struct triangulateio *object;
  int **ptriArrayI;
  int *dim0,*dim1;
  int setDim; 
  int forceResize;
  int markerDim;
  /*arguments triangulateioHandle,segments,segmarks*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin1)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArrayI = &(object->segmentmarkerlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofsegments;
  markerDim = 1;
  dim1 = &markerDim;
  /*which dimension do I need to update */
  setDim = -2; /*don't resize anything*/
  forceResize = 1; /*have to reset data arrays here*/
  rval1 = copyToTriangulateioIntegerArray1or2dim(ptriArrayI,dim0,dim1,
						 valsin1,
						 forceResize,setDim);

  return rval1;
}

/***************************************************
   get segmentlist from triangulateio data structure
   as a 2 dimensional array : Nsegments x 2
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_SEGMENTS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->segmentlist;    /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberofsegments;    /*dimensions of array to copy*/
  dim1= 2;                           /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get segmentlist from triangulateio data structure
   as a 2 dimensional array : Nsegments x 2
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_SEGMENTS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->segmentlist;    /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberofsegments;    /*dimensions of array to copy*/
  dim1= 2;                           /*dimensions of array to copy*/
  deepCopy = 1;                      /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get segmentmarkerlist from triangulateio data structure
   as a 2 dimensional array : Nsegments x 1
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_SEGMENTMARKERS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->segmentmarkerlist;    /*array to copy*/
  nd  = 1;                                 /*trianglelist is 2d array*/
  dim0= object->numberofsegments;          /*dimensions of array to copy*/
  dim1= 0;                                 /*dimensions of array to copy*/
  deepCopy = 0;                            /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get segmentmarkerlist from triangulateio data structure
   as a 2 dimensional array : Nsegments x 1
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_SEGMENTMARKERS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->segmentmarkerlist;    /*array to copy*/
  nd  = 1;                                 /*trianglelist is 2d array*/
  dim0= object->numberofsegments;          /*dimensions of array to copy*/
  dim1= 0;                                 /*dimensions of array to copy*/
  deepCopy = 1;                            /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   load hole locations into data array, assume that it is
   a 2 dimensional array : Nholes x 2
 **************************************************/
static PyObject *
triangulate_SET_HOLELIST(PyObject *self, PyObject *args){
  PyObject *address, *valsin, *rval;
  struct triangulateio *object;
  double **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int spaceDim;
  int forceResize;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->holelist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofholes;
  spaceDim = _NDIM;
  dim1 = &spaceDim;
  /*which dimension do I need to update */
  setDim = 0; /*number of attributes per triangle*/
  forceResize = 1; /*go ahead and require resizing of data arrays*/
  rval = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
					       valsin,
					       forceResize,setDim);

  /*now tell the triangulateio object that it owns memory for holelist*/
  object->holelistmemflag = 1;

  return rval;
}
/***************************************************
   get holelist from triangulateio data structure
   as a 2 dimensional array : Nholes x 2
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_HOLELIST(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->holelist;             /*array to copy*/
  nd  = 2;                                 /*trianglelist is 2d array*/
  dim0= object->numberofholes;             /*dimensions of array to copy*/
  dim1= 2;                                 /*dimensions of array to copy*/
  deepCopy = 0;                            /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}
/***************************************************
   get holelist from triangulateio data structure
   as a 2 dimensional array : Nholes x 2
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_HOLELIST_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->holelist;             /*array to copy*/
  nd  = 2;                                 /*trianglelist is 2d array*/
  dim0= object->numberofholes;             /*dimensions of array to copy*/
  dim1= 2;                                 /*dimensions of array to copy*/
  deepCopy = 1;                            /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   load regionlist constrints into data array, assume that it is
   a 2 dimensional array : Nregions x 4
   regionlist[0,0:1] are x,y coordinates, 
   regionlist[2] is regional attribute for constraint
   regionlist[3] is maximum area at constraint 
 **************************************************/
static PyObject *
triangulate_SET_REGIONLIST(PyObject *self, PyObject *args){
  PyObject *address, *valsin, *rval;
  struct triangulateio *object;
  double **ptriArray;
  int *dim0, *dim1;
  int setDim; 
  int constraintDim;
  int forceResize;
  /*arguments triangulateioHandle,attributes*/
  if(!PyArg_ParseTuple(args,(char *)"OO", 
		       &address, &valsin)){ 
    return NULL;
  }
  if(!PyCObject_Check(address)){
    PyErr_SetString(PyExc_TypeError,
      "Wrong 1st argument! CObject required (triangulateio handle).");
    return NULL;
  }    

  /*convert address to triangulateio type*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  /*pick out the array that I need to set*/
  ptriArray = &(object->regionlist);
  /*what are the correct dimensions for this list*/
  dim0 = &object->numberofregions;
  constraintDim = 4;
  dim1 = &constraintDim;
  /*which dimension do I need to update */
  setDim = 0; /*number of region/area constraints*/
  forceResize = 1; /*go ahead and require resizing of data arrays*/
  rval = copyToTriangulateioDoubleArray1or2dim(ptriArray,dim0,dim1,
					       valsin,
					       forceResize,setDim);
  /*now tell the triangulateio object that it owns memory for regionlist*/
  object->regionlistmemflag = 1;
  return rval;
}

/***************************************************
   get regionlist from triangulateio data structure
   as a 2 dimensional array : Nholes x 4
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_REGIONLIST(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->regionlist;            /*array to copy*/
  nd  = 2;                                  /*trianglelist is 2d array*/
  dim0= object->numberofregions;            /*dimensions of array to copy*/
  dim1= 4;                                  /*dimensions of array to copy*/
  deepCopy = 0;                             /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}


/***************************************************
   get regionlist from triangulateio data structure
   as a 2 dimensional array : Nholes x 4
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_REGIONLIST_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->regionlist;            /*array to copy*/
  nd  = 2;                                  /*trianglelist is 2d array*/
  dim0= object->numberofregions;            /*dimensions of array to copy*/
  dim1= 4;                                  /*dimensions of array to copy*/
  deepCopy = 1;                             /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}


/***************************************************
   get edgelist from triangulateio data structure
   as a 2 dimensional array : Nedges x 2
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_EDGELIST(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->edgelist;    /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberofedges;    /*dimensions of array to copy*/
  dim1= 2;                           /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   get edgelist from triangulateio data structure
   as a 2 dimensional array : Nedges x 2
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_EDGELIST_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->edgelist;    /*array to copy*/
  nd  = 2;                           /*trianglelist is 2d array*/
  dim0= object->numberofedges;    /*dimensions of array to copy*/
  dim1= 2;                           /*dimensions of array to copy*/
  deepCopy = 1;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}


/***************************************************
   get edgemarkerlist from triangulateio data structure
   as a 2 dimensional array : Nedges x 1
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_EDGEMARKERS(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->edgemarkerlist;    /*array to copy*/
  nd  = 2;                              /*trianglelist is 2d array*/
  dim0= object->numberofedges;          /*dimensions of array to copy*/
  dim1= 0;                              /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}


/***************************************************
   get edgemarkerlist from triangulateio data structure
   as a 2 dimensional array : Nedges x 1
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_EDGEMARKERS_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  int *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->edgemarkerlist;    /*array to copy*/
  nd  = 2;                              /*trianglelist is 2d array*/
  dim0= object->numberofedges;          /*dimensions of array to copy*/
  dim1= 0;                              /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioIntegerArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}


/***************************************************
   get normlist from triangulateio data structure
   as a 2 dimensional array : Nedges x SpaceDim
   shallow copy
 **************************************************/
static PyObject *
triangulate_GET_NORMLIST(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->normlist;   /*array to copy*/
  nd  = 2;                        /*pointslist is 2d array*/
  dim0= object->numberofedges;   /*dimensions of array to copy*/
  dim1= _NDIM;                    /*dimensions of array to copy*/
  deepCopy = 0;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}

/***************************************************
   get normlist from triangulateio data structure
   as a 2 dimensional array : Nedges x SpaceDim
   deep copy
 **************************************************/
static PyObject *
triangulate_GET_NORMLIST_COPY(PyObject *self, PyObject *args){
  PyObject *address,*rval;
  struct triangulateio *object;

  double *triArray;
  int nd,dim0,dim1,deepCopy;

  if(!PyArg_ParseTuple(args,"O", 
		       &address, &PyArray_Type))
    { 
      return NULL;
    }
  if(!PyCObject_Check(address))
    {
      PyErr_SetString(PyExc_TypeError,
		      "Wrong 1st argument! CObject required (triangulateio handle).");
      return NULL;
  }    

  /*convert address to triangulateio object*/
  object = (struct triangulateio *)PyCObject_AsVoidPtr(address);  
  triArray = object->normlist;   /*array to copy*/
  nd  = 2;                        /*pointslist is 2d array*/
  dim0= object->numberofedges;   /*dimensions of array to copy*/
  dim1= _NDIM;                    /*dimensions of array to copy*/
  deepCopy = 1;                   /*shallow or deep copy here*/
  rval = getCopyOfTriangulateioDoubleArray1or2dim(triArray,dim0,dim1,deepCopy);
  return rval;
}



/***************************************************
 document strings for methods
 **************************************************/
static char trianguleWrappers_doc[] = "provide interface to triangulateio type";
static char triangulate_NEW_doc[] = "create a new triangulateio object";
static char triangulate_APPLY_TRIANGULATE_doc[] = 
  "apply the triangulate function to input,output,and voronoi triangulateio objects";
static char triangulate_APPLY_TRIANGULATE_NO_VORONOI_doc[] = 
  "apply the triangulate function to input,output triangulateio objects omit voronoi diagram";
static char triangulate_GET_INFO_doc[] = \
  "retrieve info about triangulation: \
     numberofpoints                   \
     numberofpointattributes          \
     numberoftriangles                \
     numberoftriangleattributes       \
     numberofsegments                 \
     numberofholes                    \
     numberofregions                  \
     numberofedges";
static char triangulate_PRINT_REPORT_doc[] = 
  "just write some simple information about a triangulateio object";
static char triangulate_SET_POINTS_doc[] = 
  "set points in triangulateio: xy= Npoints x 2";
static char triangulate_SET_POINTS_AND_MARKERS_doc[] = 
  "set points and markers in triangulateio: xy= Npoints x 2, markers= Npoints ";
static char triangulate_SET_POINTMARKERS_doc[] = 
  "set pointmarkers in triangulateio:markers= Npoints ";
static char triangulate_GET_POINTS_COPY_doc[] = 
  "retrieve points in triangulateio: xy= Npoints x 2 deep copy ";
static char triangulate_GET_POINTS_doc[] = 
  "retrieve points in triangulateio: xy= Npoints x 2 ";
static char triangulate_GET_POINTMARKERS_doc[] = 
  "retrieve pointmarkers in triangulateio: xy= Npoints x 1 ";
static char triangulate_GET_POINTMARKERS_COPY_doc[] = 
  "retrieve copy of pointmarkers in triangulateio: xy= Npoints x 1 ";
static char triangulate_SET_POINT_ATTRIBS_doc[] = 
  "set attributes of points in triangulateio: a= Npoints x Nattrib ";
static char triangulate_GET_POINT_ATTRIBS_doc[] = 
  "get attributes of triangles in triangulateio: a= Nelems x Nattrib ";
static char triangulate_GET_POINT_ATTRIBS_COPY_doc[] = 
  "get copy attributes of triangles in triangulateio: a= Nelems x Nattrib ";
static char triangulate_SET_TRIANGLES_doc[] = 
  "set triangles in triangulateio: xy= Nelem x Ncorner (3 or 6)";
static char triangulate_GET_TRIANGLES_doc[] = 
  "get triangles in triangulateio: tri= Nelem x Ncorner (3 or 6)";
static char triangulate_GET_TRIANGLES_COPY_doc[] = 
  "get copy of triangles in triangulateio: tri= Nelem x Ncorner (3 or 6)";
static char triangulate_SET_TRIANGLE_ATTRIBS_doc[] = 
  "set attributes of triangles in triangulateio: a= Nelems x Nattrib ";
static char triangulate_GET_TRIANGLE_ATTRIBS_doc[] = 
  "get attributes of triangles in triangulateio: a= Nelems x Nattrib ";
static char triangulate_GET_TRIANGLE_ATTRIBS_COPY_doc[] = 
  "get copy attributes of triangles in triangulateio: a= Nelems x Nattrib ";
static char triangulate_SET_TRIANGLE_AREAS_doc[] = 
  "set area contraints of triangles in triangulateio: a= Nelems x 1 ";
static char triangulate_GET_NEIGHBORS_doc[] = 
  "get triangles in triangulateio: tri= Nelem x 3";
static char triangulate_GET_NEIGHBORS_COPY_doc[] = 
  "get copy of triangles in triangulateio: tri= Nelem x 3";
static char triangulate_SET_SEGMENTS_doc[] = 
  "set boundary in triangulateio: xy= Nsegments x 2";
static char triangulate_SET_SEGMENTS_AND_MARKERS_doc[] = 
  "set boundary in triangulateio: xy= Nsegments x 2, markers= Nsegments ";
static char triangulate_SET_SEGMENTMARKERS_doc[] = 
  "set boundary markers in triangulateio:  markers= Nsegments ";
static char triangulate_GET_SEGMENTS_doc[] = 
  "get boundary segments in triangulateio: xy= Nsegments x 2";
static char triangulate_GET_SEGMENTS_COPY_doc[] = 
  "get copy of boundary segments in triangulateio: xy= Nsegments x 2";
static char triangulate_GET_SEGMENTMARKERS_doc[] = 
  "get boundary segmentmarkers in triangulateio: xy= Nsegments x 1";
static char triangulate_GET_SEGMENTMARKERS_COPY_doc[] = 
  "get copy of boundary segmentmarkers in triangulateio: xy= Nsegments x 1";
static char triangulate_SET_HOLELIST_doc[] = 
  "set holes in triangulateio: xy= Nholes x 2 ";
static char triangulate_GET_HOLELIST_doc[] = 
  "get holes in triangulateio: xy= Nholes x 2 ";
static char triangulate_GET_HOLELIST_COPY_doc[] = 
  "get copy of holes in triangulateio: xy= Nholes x 2 ";
static char triangulate_SET_REGIONLIST_doc[] = 
  "set regionlist/area constraints in triangulateio: xy= Nregions x 4 ";
static char triangulate_GET_REGIONLIST_doc[] = 
  "get regionlist/area constraints in triangulateio: xy= Nregions x 4 ";
static char triangulate_GET_REGIONLIST_COPY_doc[] = 
  "get copy of regionlist/area constraints in triangulateio: xy= Nregions x 4 ";
static char triangulate_GET_EDGELIST_doc[] = 
  "get edgelist in triangulateio: xy= Nedges x 2 ";
static char triangulate_GET_EDGELIST_COPY_doc[] = 
  "get copy of edgelist in triangulateio: xy= Nedges x 2 ";
static char triangulate_GET_EDGEMARKERS_doc[] = 
  "get edgemarkers in triangulateio: xy= Nedges x 2 ";
static char triangulate_GET_EDGEMARKERS_COPY_doc[] = 
  "get copy of edgemarkers in triangulateio: xy= Nedges x 2 ";
static char triangulate_GET_NORMLIST_doc[] = 
  "get normlist in triangulateio: xy= Nedges x 2 ";
static char triangulate_GET_NORMLIST_COPY_doc[] = 
  "get copy of normlist in triangulateio: xy= Nedges x 2 ";

/***************************************************
 method table
 **************************************************/
static PyMethodDef triangleWrappersMethods[] = {
  {"new",                /*name when called from python*/
   triangulate_NEW,      /*C function name */
   METH_VARARGS,         /*ordinary (no keyword) arguments */
   triangulate_NEW_doc}, /*doc string for method*/
  {"applyTriangulate",   /*name when called from python*/
   triangulate_APPLY_TRIANGULATE,      /*C function name */
   METH_VARARGS,                       /*ordinary (no keyword) arguments */
   triangulate_APPLY_TRIANGULATE_doc}, /*doc string for method*/
  {"applyTriangulateNoVoronoi",                   /*name when called from python*/
   triangulate_APPLY_TRIANGULATE_NO_VORONOI,      /*C function name */
   METH_VARARGS,                                  /*ordinary (no keyword) arguments */
   triangulate_APPLY_TRIANGULATE_NO_VORONOI_doc}, /*doc string for method*/
  {"getInfo",                   /*name when called from python*/
   triangulate_GET_INFO,        /*C function name */
   METH_VARARGS,                /*ordinary (no keyword) arguments */
   triangulate_GET_INFO_doc},       /*doc string for method*/
  {"printReport",                   /*name when called from python*/
   triangulate_PRINT_REPORT,        /*C function name */
   METH_VARARGS,                    /*ordinary (no keyword) arguments */
   triangulate_PRINT_REPORT_doc},   /*doc string for method*/
  {"setPoints",                 /*name when called from python*/
   triangulate_SET_POINTS,      /*C function name */
   METH_VARARGS,                /*ordinary (no keyword) arguments */
   triangulate_SET_POINTS_doc}, /*doc string for method*/
  {"setPointsAndMarkers",                   /*name when called from python*/
   triangulate_SET_POINTS_AND_MARKERS,      /*C function name */
   METH_VARARGS,                            /*ordinary (no keyword) arguments */
   triangulate_SET_POINTS_AND_MARKERS_doc}, /*doc string for method*/
  {"setPointMarkers",                 /*name when called from python*/
   triangulate_SET_POINTMARKERS,      /*C function name */
   METH_VARARGS,                      /*ordinary (no keyword) arguments */
   triangulate_SET_POINTMARKERS_doc}, /*doc string for method*/
  {"getPoints",                 /*name when called from python*/
   triangulate_GET_POINTS,      /*C function name */
   METH_VARARGS,                /*ordinary (no keyword) arguments */
   triangulate_GET_POINTS_doc},      /*doc string for method*/
  {"getPointsCopy",                  /*name when called from python*/
   triangulate_GET_POINTS_COPY,      /*C function name */
   METH_VARARGS,                     /*ordinary (no keyword) arguments */
   triangulate_GET_POINTS_COPY_doc}, /*doc string for method*/
  {"getPointMarkers",                 /*name when called from python*/
   triangulate_GET_POINTMARKERS,      /*C function name */
   METH_VARARGS,                /*ordinary (no keyword) arguments */
   triangulate_GET_POINTMARKERS_doc},      /*doc string for method*/
  {"getPointMarkersCopy",                 /*name when called from python*/
   triangulate_GET_POINTMARKERS_COPY,      /*C function name */
   METH_VARARGS,                /*ordinary (no keyword) arguments */
   triangulate_GET_POINTMARKERS_COPY_doc},      /*doc string for method*/
  {"setPointAttributes",                /*name when called from python*/
   triangulate_SET_POINT_ATTRIBS,       /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_SET_POINT_ATTRIBS_doc},  /*doc string for method*/
  {"getPointAttributes",                /*name when called from python*/
   triangulate_GET_POINT_ATTRIBS,       /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_POINT_ATTRIBS_doc},  /*doc string for method*/
  {"getPointAttributesCopy",                 /*name when called from python*/
   triangulate_GET_POINT_ATTRIBS_COPY,       /*C function name */
   METH_VARARGS,                             /*ordinary (no keyword) arguments */
   triangulate_GET_POINT_ATTRIBS_COPY_doc},  /*doc string for method*/
  {"setTriangles",                 /*name when called from python*/
   triangulate_SET_TRIANGLES,      /*C function name */
   METH_VARARGS,                  /*ordinary (no keyword) arguments */
   triangulate_SET_TRIANGLES_doc}, /*doc string for method*/
  {"getTriangles",                 /*name when called from python*/
   triangulate_GET_TRIANGLES,      /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_TRIANGLES_doc}, /*doc string for method*/
  {"getTrianglesCopy",                  /*name when called from python*/
   triangulate_GET_TRIANGLES_COPY,      /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_TRIANGLES_COPY_doc}, /*doc string for method*/
  {"setTriangleAttributes",                /*name when called from python*/
   triangulate_SET_TRIANGLE_ATTRIBS,       /*C function name */
   METH_VARARGS,                           /*ordinary (no keyword) arguments */
   triangulate_SET_TRIANGLE_ATTRIBS_doc},  /*doc string for method*/
  {"getTriangleAttributes",                /*name when called from python*/
   triangulate_GET_TRIANGLE_ATTRIBS,       /*C function name */
   METH_VARARGS,                           /*ordinary (no keyword) arguments */
   triangulate_GET_TRIANGLE_ATTRIBS_doc}, /*doc string for method*/
  {"getTriangleAttributesCopy",                 /*name when called from python*/
   triangulate_GET_TRIANGLE_ATTRIBS_COPY,       /*C function name */
   METH_VARARGS,                                /*ordinary (no keyword) arguments */
   triangulate_GET_TRIANGLE_ATTRIBS_COPY_doc}, /*doc string for method*/
  {"setTriangleAreas",                    /*name when called from python*/
   triangulate_SET_TRIANGLE_AREAS,        /*C function name */
   METH_VARARGS,                          /*ordinary (no keyword) arguments */
   triangulate_SET_TRIANGLE_AREAS_doc},   /*doc string for method*/
  {"getNeighbors",                 /*name when called from python*/
   triangulate_GET_NEIGHBORS,      /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_NEIGHBORS_doc}, /*doc string for method*/
  {"getNeighborsCopy",                  /*name when called from python*/
   triangulate_GET_NEIGHBORS_COPY,      /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_NEIGHBORS_COPY_doc}, /*doc string for method*/
  {"setSegments",                  /*name when called from python*/
   triangulate_SET_SEGMENTS,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_SET_SEGMENTS_doc},  /*doc string for method*/
  {"setSegmentsAndMarkers",                   /*name when called from python*/
   triangulate_SET_SEGMENTS_AND_MARKERS,      /*C function name */
   METH_VARARGS,                              /*ordinary (no keyword) arguments */
   triangulate_SET_SEGMENTS_AND_MARKERS_doc}, /*doc string for method*/
  {"setSegmentMarkers",                  /*name when called from python*/
   triangulate_SET_SEGMENTMARKERS,       /*C function name */
   METH_VARARGS,                         /*ordinary (no keyword) arguments */
   triangulate_SET_SEGMENTMARKERS_doc},  /*doc string for method*/
  {"getSegments",                  /*name when called from python*/
   triangulate_GET_SEGMENTS,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_SEGMENTS_doc},  /*doc string for method*/
  {"getSegmentsCopy",                   /*name when called from python*/
   triangulate_GET_SEGMENTS_COPY,       /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_SEGMENTS_COPY_doc},  /*doc string for method*/
  {"getSegmentMarkers",                  /*name when called from python*/
   triangulate_GET_SEGMENTMARKERS,       /*C function name */
   METH_VARARGS,                         /*ordinary (no keyword) arguments */
   triangulate_GET_SEGMENTMARKERS_doc}, /*doc string for method*/
  {"getSegmentMarkersCopy",                   /*name when called from python*/
   triangulate_GET_SEGMENTMARKERS_COPY,       /*C function name */
   METH_VARARGS,                              /*ordinary (no keyword) arguments */
   triangulate_GET_SEGMENTMARKERS_COPY_doc}, /*doc string for method*/
  {"setHoles",                     /*name when called from python*/
   triangulate_SET_HOLELIST,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_SET_HOLELIST_doc},  /*doc string for method*/
  {"getHoles",                     /*name when called from python*/
   triangulate_GET_HOLELIST,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_HOLELIST_doc},  /*doc string for method*/
  {"getHolesCopy",                     /*name when called from python*/
   triangulate_GET_HOLELIST_COPY,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_HOLELIST_COPY_doc},  /*doc string for method*/
  {"setRegions",                     /*name when called from python*/
   triangulate_SET_REGIONLIST,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_SET_REGIONLIST_doc},  /*doc string for method*/
  {"getRegions",                     /*name when called from python*/
   triangulate_GET_REGIONLIST,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_REGIONLIST_doc},  /*doc string for method*/
  {"getRegionsCopy",                     /*name when called from python*/
   triangulate_GET_REGIONLIST_COPY,       /*C function name */
   METH_VARARGS,                          /*ordinary (no keyword) arguments */
   triangulate_GET_REGIONLIST_COPY_doc},  /*doc string for method*/
  {"getEdges",                     /*name when called from python*/
   triangulate_GET_EDGELIST,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_EDGELIST_doc},  /*doc string for method*/
  {"getEdgesCopy",                      /*name when called from python*/
   triangulate_GET_EDGELIST_COPY,       /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_EDGELIST_COPY_doc},  /*doc string for method*/
  {"getEdgeMarkers",                    /*name when called from python*/
   triangulate_GET_EDGEMARKERS,         /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_EDGEMARKERS_doc},    /*doc string for method*/
  {"getEdgeMarkersCopy",                   /*name when called from python*/
   triangulate_GET_EDGEMARKERS_COPY,       /*C function name */
   METH_VARARGS,                           /*ordinary (no keyword) arguments */
   triangulate_GET_EDGEMARKERS_COPY_doc},  /*doc string for method*/
   {"getVoronoiNormals",           /*name when called from python*/
   triangulate_GET_NORMLIST,       /*C function name */
   METH_VARARGS,                   /*ordinary (no keyword) arguments */
   triangulate_GET_NORMLIST_doc},  /*doc string for method*/
   {"getVoronoiNormalsCopy",            /*name when called from python*/
   triangulate_GET_NORMLIST_COPY,       /*C function name */
   METH_VARARGS,                        /*ordinary (no keyword) arguments */
   triangulate_GET_NORMLIST_COPY_doc},  /*doc string for method*/
  {NULL,NULL}                             /*required ending for method table*/
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC inittriangleWrappers(void)
{
  PyObject *m,*d;
  m = Py_InitModule3("triangleWrappers", triangleWrappersMethods,trianguleWrappers_doc);
  d = PyModule_GetDict(m);
  import_array();
}
/** @} */
