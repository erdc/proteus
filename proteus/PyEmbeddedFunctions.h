#ifndef PYEMBEDDEDFUNCTIONS_H
#define PYEMBEDDEDFUNCTIONS_H

#include "Python.h"

//This function accepts a string which gets passed to the logEvent Python function in Profiling
int logEvent(char* logString,int logLevel)

//logString is the desired string: usually set with sprintf() to combine characters with numbers
//logLevel means the same as with the python function
{
  if(!Py_IsInitialized())
    Py_Initialize();
  PyObject *pName=NULL, *pModule=NULL, *pFunc=NULL;
  PyObject *pArgs=NULL, *pValue=NULL;

  pName = PyUnicode_FromString("proteus.Profiling");
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (pModule != NULL)
  {
    pFunc = PyObject_GetAttrString(pModule,"logEvent");
    pArgs = PyTuple_New(2);

/*
  //This is how the embedding is done in the official Python example...but this leads to a memory leak for some reason
    pValue = PyUnicode_FromString(logString);
    // pValue reference stolen here:
    PyTuple_SetItem(pArgs, 0, pValue);

    pValue = PyInt_FromLong(logLevel);
    // pValue reference stolen here:
    PyTuple_SetItem(pArgs, 1, pValue);
*/
    PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(logString));
    PyTuple_SetItem(pArgs, 1, PyInt_FromLong(logLevel));

    PyObject_CallObject(pFunc,pArgs);
    Py_DECREF(pFunc);
    Py_DECREF(pArgs);
    if(pValue !=NULL)
    {
      Py_DECREF(pValue);
      Py_XDECREF(pValue);
    }
  }
  else
  {
    PyErr_Print();
    fprintf(stderr,"Failed to load \"%s\"\n", "proteus.Profiling"); 
    return 1;
  }
  Py_DECREF(pModule);
  return 0;
}

#endif
