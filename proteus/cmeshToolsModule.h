#ifndef CMESHTOOLS_MODULE_H
#define CMESHTOOLS_MODULE_H

#include "Python.h"
#include "numpy/arrayobject.h"
#include "mesh.h"

typedef struct
{
  PyObject_HEAD
  Mesh mesh;
} CMesh;

#define MESH(p) ((CMesh*)p)->mesh





#endif
