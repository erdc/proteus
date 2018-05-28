#ifndef PYEMBEDDEDFUNCTIONS_H
#define PYEMBEDDEDFUNCTIONS_H

#include "Python.h"

//I want this function to accept a string which gets passed to the logEvent Python function in Profiling
int logEvent(char logString[],int logLevel)
{
  if(!Py_IsInitialized())
    Py_Initialize();
  //PyRun_SimpleString("print 'Test Embedding'\n");
  PyObject *pName, *pModule, *pFunc;
  PyObject *pArgs, *pValue;

  pName = PyString_FromString("proteus.Profiling");
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);
  if (pModule != NULL)
  {
    pFunc = PyObject_GetAttrString(pModule,"logEvent");
    pArgs = PyTuple_New(2);

    pValue = PyString_FromString(logString);
    /* pValue reference stolen here: */
    PyTuple_SetItem(pArgs, 0, pValue);

    pValue = PyInt_FromLong(logLevel);
    /* pValue reference stolen here: */

    PyTuple_SetItem(pArgs, 1, pValue);
    PyObject_CallObject(pFunc,pArgs);
    Py_DECREF(pFunc);
    Py_DECREF(pArgs);
    if(pValue !=NULL)
    {
      Py_DECREF(pValue);
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
