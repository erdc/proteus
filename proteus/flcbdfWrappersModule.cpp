/** \file flcbdfWrappersModule.cpp
    \ingroup flcbdfWrappers
    @{
*/
#define FLCBDF_WRAPPERS_MODULE
#include "flcbdfWrappersModule.h"
#include <algorithm>
#include "meshio.h"

using namespace proteus;

#define DDATA(p) ((double *) (((PyArrayObject *)p)->data))
#define SHAPE(p) ((PyArrayObject *)p)->dimensions
#define ND(p) ((PyArrayObject *)p)->nd
#define IDATA(p) ((int *) (((PyArrayObject *)p)->data))
#define SMP(p) ((SparseMatrix*)p)
#define MESH(p) ((CMesh*)p)->mesh


extern "C"
{
  static int
  ensure_comm()
  {
    return 1;
    // if (PROTEUS_COMM_WORLD == MPI_COMM_NULL) {
    //   PyErr_SetString(PyExc_RuntimeError, "flcbdfWrappersModule is not initialized!");
    //   return 0;
    // }
    // return 1;
  }

  static PyObject* flcbdfWrappersGlobalSum(PyObject* self, PyObject* args)
  {
    using namespace std;
    double value,value_new;

    if (!ensure_comm()) {
      return NULL;
    }

    if (!PyArg_ParseTuple(args,
                          "d",
                          &value))
      return NULL;
    MPI_Allreduce(&value,&value_new,1,MPI_DOUBLE,MPI_SUM,PROTEUS_COMM_WORLD);
    return Py_BuildValue("d",value_new);
  }

  static PyObject* flcbdfWrappersGlobalMax(PyObject* self, PyObject* args)
  {
    using namespace std;
    double value,value_new;

    if (!ensure_comm()) {
      return NULL;
    }

    if (!PyArg_ParseTuple(args,
                          "d",
                          &value))
      return NULL;

    MPI_Allreduce(&value,&value_new,1,MPI_DOUBLE,MPI_MAX,PROTEUS_COMM_WORLD);
    return Py_BuildValue("d",value_new);
  }
  
  static PyObject* flcbdfWrappersGlobalMin(PyObject* self, PyObject* args)
  {
    using namespace std;
    double value,value_new;

    if (!ensure_comm()) {
      return NULL;
    }

    if (!PyArg_ParseTuple(args,
                          "d",
                          &value))
      return NULL;
    MPI_Allreduce(&value,&value_new,1,MPI_DOUBLE,MPI_MIN,PROTEUS_COMM_WORLD);
    return Py_BuildValue("d",value_new);
  }

  static PyObject* flcbdfWrappersPartitionElements(PyObject* self,
                                                   PyObject* args)
  {
    using namespace std;
    int nLayersOfOverlap;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *nodeOffsets_subdomain_owned,
      *nodeNumbering_subdomain2global,
      *elementBoundaryOffsets_subdomain_owned,
      *elementBoundaryNumbering_subdomain2global,
      *edgeOffsets_subdomain_owned,
      *edgeNumbering_subdomain2global;
    if (!PyArg_ParseTuple(args,
                          "iOO",
                          &nLayersOfOverlap,
                          &cmesh,
                          &subdomain_cmesh))
      return NULL;
    MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
    PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;
    int ierr,size,rank;

    if (!ensure_comm()) {
      return NULL;
    }

    ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
    partitionElements(PROTEUS_COMM_WORLD, MESH(cmesh),nLayersOfOverlap);

    int dims[1];
    //build handles to python arrays
    dims[0] = size+1;
    elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).elementOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElements_global;
    elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                dims,
                                                                PyArray_INT,
                                                                (char*)MESH(cmesh).elementNumbering_subdomain2global);
    dims[0] = size+1;
    nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).nodeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nNodes_global;
    nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).nodeNumbering_subdomain2global);
    dims[0] = size+1;
    elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                     dims,
                                                                     PyArray_INT,
                                                                     (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
    elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                        dims,
                                                                        PyArray_INT,
                                                                        (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
    dims[0] = size+1;
    edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).edgeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nEdges_global;
    edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).edgeNumbering_subdomain2global);
    return Py_BuildValue("OOOOOOOO",
                         elementOffsets_subdomain_owned,
                         elementNumbering_subdomain2global,
                         nodeOffsets_subdomain_owned,
                         nodeNumbering_subdomain2global,
                         elementBoundaryOffsets_subdomain_owned,
                         elementBoundaryNumbering_subdomain2global,
                         edgeOffsets_subdomain_owned,
                         edgeNumbering_subdomain2global);
  }

  static PyObject* flcbdfWrappersPartitionNodes(PyObject* self,
                                                PyObject* args)
  {
    using namespace std;
    int nLayersOfOverlap;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *nodeOffsets_subdomain_owned,
      *nodeNumbering_subdomain2global,
      *elementBoundaryOffsets_subdomain_owned,
      *elementBoundaryNumbering_subdomain2global,
      *edgeOffsets_subdomain_owned,
      *edgeNumbering_subdomain2global;
    if (!PyArg_ParseTuple(args,
                          "iOO",
                          &nLayersOfOverlap,
                          &cmesh,
                          &subdomain_cmesh))
      return NULL;

    if (!ensure_comm()) {
      return NULL;
    }

    MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
    PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;
    int ierr,size,rank;
    ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
    partitionNodes(PROTEUS_COMM_WORLD, MESH(cmesh),nLayersOfOverlap);

    int dims[1];
    //build handles to python arrays
    dims[0] = size+1;
    elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).elementOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElements_global;
    elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                dims,
                                                                PyArray_INT,
                                                                (char*)MESH(cmesh).elementNumbering_subdomain2global);
    dims[0] = size+1;
    nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).nodeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nNodes_global;
    nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).nodeNumbering_subdomain2global);
    dims[0] = size+1;
    elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                     dims,
                                                                     PyArray_INT,
                                                                     (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
    elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                        dims,
                                                                        PyArray_INT,
                                                                        (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
    dims[0] = size+1;
    edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).edgeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nEdges_global;
    edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).edgeNumbering_subdomain2global);

    return Py_BuildValue("OOOOOOOO",
                         elementOffsets_subdomain_owned,
                         elementNumbering_subdomain2global,
                         nodeOffsets_subdomain_owned,
                         nodeNumbering_subdomain2global,
                         elementBoundaryOffsets_subdomain_owned,
                         elementBoundaryNumbering_subdomain2global,
                         edgeOffsets_subdomain_owned,
                         edgeNumbering_subdomain2global);
  }

  static PyObject* flcbdfWrappersConvertPUMIPartitionToPython(PyObject* self,
                                                              PyObject* args)
  {
    using namespace std;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *elementNumbering_global2original,
      *nodeOffsets_subdomain_owned,
      *nodeNumbering_subdomain2global,
      *nodeNumbering_global2original,
      *elementBoundaryOffsets_subdomain_owned,
      *elementBoundaryNumbering_subdomain2global,
      *elementBoundaryNumbering_global2original,
      *edgeOffsets_subdomain_owned,
      *edgeNumbering_subdomain2global,
      *edgeNumbering_global2original;
    if (!PyArg_ParseTuple(args,
                          "OO",
                          &cmesh,
                          &subdomain_cmesh))
      return NULL;
    MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
    PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;
    int ierr,size,rank;
    ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);

    int dims[1];
    //build handles to python arrays
    dims[0] = size+1;
    elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).elementOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElements_global;
    elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                dims,
                                                                PyArray_INT,
                                                                (char*)MESH(cmesh).elementNumbering_subdomain2global);
    dims[0] = size+1;
    nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).nodeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nNodes_global;
    nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).nodeNumbering_subdomain2global);
    dims[0] = size+1;
    elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                     dims,
                                                                     PyArray_INT,
                                                                     (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
    elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                        dims,
                                                                        PyArray_INT,
                                                                        (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
    dims[0] = size+1;
    edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).edgeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nEdges_global;
    edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).edgeNumbering_subdomain2global);
    return Py_BuildValue("OOOOOOOO",
                         elementOffsets_subdomain_owned,
                         elementNumbering_subdomain2global,
                         nodeOffsets_subdomain_owned,
                         nodeNumbering_subdomain2global,
                         elementBoundaryOffsets_subdomain_owned,
                         elementBoundaryNumbering_subdomain2global,
                         edgeOffsets_subdomain_owned,
                         edgeNumbering_subdomain2global);

  }

  static PyObject* flcbdfWrappersPartitionNodesFromTetgenFiles(PyObject* self,
                                                               PyObject* args)
  {
    using namespace std;
    int nLayersOfOverlap, indexBase;
    char* filebase;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *nodeOffsets_subdomain_owned,
      *nodeNumbering_subdomain2global,
      *elementBoundaryOffsets_subdomain_owned,
      *elementBoundaryNumbering_subdomain2global,
      *edgeOffsets_subdomain_owned,
      *edgeNumbering_subdomain2global;
    if (!PyArg_ParseTuple(args,
                          "siiOO",
                          &filebase,
                          &indexBase,
                          &nLayersOfOverlap,
                          &cmesh,
                          &subdomain_cmesh))
      return NULL;
    MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
    PETSC_COMM_WORLD = PROTEUS_COMM_WORLD;

    if (!ensure_comm()) {
      return NULL;
    }

    int ierr,size,rank;
    ierr = MPI_Comm_size(PROTEUS_COMM_WORLD,&size);
    ierr = MPI_Comm_rank(PROTEUS_COMM_WORLD,&rank);
    partitionNodesFromTetgenFiles(PROTEUS_COMM_WORLD, filebase,indexBase,MESH(cmesh),nLayersOfOverlap);

    int dims[1];
    //build handles to python arrays
    dims[0] = size+1;
    elementOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).elementOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElements_global;
    elementNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                dims,
                                                                PyArray_INT,
                                                                (char*)MESH(cmesh).elementNumbering_subdomain2global);
    dims[0] = size+1;
    nodeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).nodeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nNodes_global;
    nodeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).nodeNumbering_subdomain2global);
    dims[0] = size+1;
    elementBoundaryOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                                     dims,
                                                                     PyArray_INT,
                                                                     (char*)MESH(cmesh).elementBoundaryOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nElementBoundaries_global;
    elementBoundaryNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                                        dims,
                                                                        PyArray_INT,
                                                                        (char*)MESH(cmesh).elementBoundaryNumbering_subdomain2global);
    dims[0] = size+1;
    edgeOffsets_subdomain_owned = PyArray_FromDimsAndData(1,
                                                          dims,
                                                          PyArray_INT,
                                                          (char*)MESH(cmesh).edgeOffsets_subdomain_owned);

    dims[0] = MESH(cmesh).subdomainp->nEdges_global;
    edgeNumbering_subdomain2global = PyArray_FromDimsAndData(1,
                                                             dims,
                                                             PyArray_INT,
                                                             (char*)MESH(cmesh).edgeNumbering_subdomain2global);
    return Py_BuildValue("OOOOOOOO",
                         elementOffsets_subdomain_owned,
                         elementNumbering_subdomain2global,
                         nodeOffsets_subdomain_owned,
                         nodeNumbering_subdomain2global,
                         elementBoundaryOffsets_subdomain_owned,
                         elementBoundaryNumbering_subdomain2global,
                         edgeOffsets_subdomain_owned,
                         edgeNumbering_subdomain2global);
  }

  static PyObject* flcbdfWrappersBuildQuadraticLocal2GlobalMappings(PyObject* self,
                                                                    PyObject* args)
  {
    using namespace std;
    int nSpace;
    int nDOF_all_processes=0, nDOF_subdomain=0, max_dof_neighbors=0;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *nodeOffsets_subdomain_owned,
      *elementBoundaryOffsets_subdomain_owned,
      *edgeOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *nodeNumbering_subdomain2global,
      *elementBoundaryNumbering_subdomain2global,
      *edgeNumbering_subdomain2global,
      *quadratic_dof_offsets_subdomain_owned,
      *quadraticNumbering_subdomain2global,
      *quadratic_subdomain_l2g,
      *quadratic_lagrangeNodes;
    if (!PyArg_ParseTuple(args,
                          "iOOOOOOOOOOOOOO",
                          &nSpace,
                          &cmesh,
                          &subdomain_cmesh,
                          &elementOffsets_subdomain_owned,
                          &nodeOffsets_subdomain_owned,
                          &elementBoundaryOffsets_subdomain_owned,
                          &edgeOffsets_subdomain_owned,
                          &elementNumbering_subdomain2global,
                          &nodeNumbering_subdomain2global,
                          &elementBoundaryNumbering_subdomain2global,
                          &edgeNumbering_subdomain2global,
                          &quadratic_dof_offsets_subdomain_owned,
                          &quadratic_subdomain_l2g,
                          &quadraticNumbering_subdomain2global,
                          &quadratic_lagrangeNodes))
      return NULL;
    //MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
    if (nSpace == 1)
      {
        buildQuadraticSubdomain2GlobalMappings_1d(PROTEUS_COMM_WORLD,
                                                  MESH(cmesh),
                                                  IDATA(elementOffsets_subdomain_owned),
                                                  IDATA(nodeOffsets_subdomain_owned),
                                                  IDATA(elementNumbering_subdomain2global),
                                                  IDATA(nodeNumbering_subdomain2global),
                                                  nDOF_all_processes,
                                                  nDOF_subdomain,
                                                  max_dof_neighbors,
                                                  IDATA(quadratic_dof_offsets_subdomain_owned),
                                                  IDATA(quadratic_subdomain_l2g),
                                                  IDATA(quadraticNumbering_subdomain2global),
                                                  DDATA(quadratic_lagrangeNodes));
      }
    else if (nSpace == 2)
      {
        buildQuadraticSubdomain2GlobalMappings_2d(PROTEUS_COMM_WORLD,
                                                  MESH(cmesh),
                                                  IDATA(elementBoundaryOffsets_subdomain_owned),
                                                  IDATA(nodeOffsets_subdomain_owned),
                                                  IDATA(elementBoundaryNumbering_subdomain2global),
                                                  IDATA(nodeNumbering_subdomain2global),
                                                  nDOF_all_processes,
                                                  nDOF_subdomain,
                                                  max_dof_neighbors,
                                                  IDATA(quadratic_dof_offsets_subdomain_owned),
                                                  IDATA(quadratic_subdomain_l2g),
                                                  IDATA(quadraticNumbering_subdomain2global),
                                                  DDATA(quadratic_lagrangeNodes));

      }
    else
      {
        buildQuadraticSubdomain2GlobalMappings_3d(PROTEUS_COMM_WORLD,
                                                  MESH(cmesh),
                                                  IDATA(edgeOffsets_subdomain_owned),
                                                  IDATA(nodeOffsets_subdomain_owned),
                                                  IDATA(edgeNumbering_subdomain2global),
                                                  IDATA(nodeNumbering_subdomain2global),
                                                  nDOF_all_processes,
                                                  nDOF_subdomain,
                                                  max_dof_neighbors,
                                                  IDATA(quadratic_dof_offsets_subdomain_owned),
                                                  IDATA(quadratic_subdomain_l2g),
                                                  IDATA(quadraticNumbering_subdomain2global),
                                                  DDATA(quadratic_lagrangeNodes));
      }


    return Py_BuildValue("iii",
                         nDOF_all_processes,
                         nDOF_subdomain,
                         max_dof_neighbors);

  }

  static PyObject* flcbdfWrappersBuildQuadraticCubeLocal2GlobalMappings(PyObject* self,
                                                                        PyObject* args)
  {
    using namespace std;
    int nSpace;
    int nDOF_all_processes=0, nDOF_subdomain=0, max_dof_neighbors=0;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *nodeOffsets_subdomain_owned,
      *elementBoundaryOffsets_subdomain_owned,
      *edgeOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *nodeNumbering_subdomain2global,
      *elementBoundaryNumbering_subdomain2global,
      *edgeNumbering_subdomain2global,
      *quadratic_dof_offsets_subdomain_owned,
      *quadraticNumbering_subdomain2global,
      *quadratic_subdomain_l2g,
      *quadratic_lagrangeNodes;
    if (!PyArg_ParseTuple(args,
                          "iOOOOOOOOOOOOOO",
                          &nSpace,
                          &cmesh,
                          &subdomain_cmesh,
                          &elementOffsets_subdomain_owned,
                          &nodeOffsets_subdomain_owned,
                          &elementBoundaryOffsets_subdomain_owned,
                          &edgeOffsets_subdomain_owned,
                          &elementNumbering_subdomain2global,
                          &nodeNumbering_subdomain2global,
                          &elementBoundaryNumbering_subdomain2global,
                          &edgeNumbering_subdomain2global,
                          &quadratic_dof_offsets_subdomain_owned,
                          &quadratic_subdomain_l2g,
                          &quadraticNumbering_subdomain2global,
                          &quadratic_lagrangeNodes))
      return NULL;
    //MESH(cmesh).subdomainp=&MESH(subdomain_cmesh);
    if (nSpace == 1)
      {
        /*buildQuadraticCubeSubdomain2GlobalMappings_1d(MESH(cmesh),
          IDATA(elementOffsets_subdomain_owned),
          IDATA(nodeOffsets_subdomain_owned),
          IDATA(elementNumbering_subdomain2global),
          IDATA(nodeNumbering_subdomain2global),
          nDOF_all_processes,
          nDOF_subdomain,
          max_dof_neighbors,
          IDATA(quadratic_dof_offsets_subdomain_owned),
          IDATA(quadratic_subdomain_l2g),
          IDATA(quadraticNumbering_subdomain2global),
          DDATA(quadratic_lagrangeNodes));*/
        std::cout<<"buildQuadraticCubeSubdomain2GlobalMappings_1d not implemented!!"<<std::endl;
      }
    else if (nSpace == 2)
      {
        /* buildQuadraticCubeSubdomain2GlobalMappings_2d(MESH(cmesh),
           IDATA(elementBoundaryOffsets_subdomain_owned),
           IDATA(nodeOffsets_subdomain_owned),
           IDATA(elementBoundaryNumbering_subdomain2global),
           IDATA(nodeNumbering_subdomain2global),
           nDOF_all_processes,
           nDOF_subdomain,
           max_dof_neighbors,
           IDATA(quadratic_dof_offsets_subdomain_owned),
           IDATA(quadratic_subdomain_l2g),
           IDATA(quadraticNumbering_subdomain2global),
           DDATA(quadratic_lagrangeNodes));*/
        std::cout<<"buildQuadraticCubeSubdomain2GlobalMappings_2d not implemented!!"<<std::endl;
      }
    else
      {
        buildQuadraticCubeSubdomain2GlobalMappings_3d(PROTEUS_COMM_WORLD,
                                                      MESH(cmesh),
                                                      IDATA(edgeOffsets_subdomain_owned),
                                                      IDATA(nodeOffsets_subdomain_owned),
                                                      IDATA(edgeNumbering_subdomain2global),
                                                      IDATA(nodeNumbering_subdomain2global),
                                                      nDOF_all_processes,
                                                      nDOF_subdomain,
                                                      max_dof_neighbors,
                                                      IDATA(quadratic_dof_offsets_subdomain_owned),
                                                      IDATA(quadratic_subdomain_l2g),
                                                      IDATA(quadraticNumbering_subdomain2global),
                                                      DDATA(quadratic_lagrangeNodes));
      }


    return Py_BuildValue("iii",
                         nDOF_all_processes,
                         nDOF_subdomain,
                         max_dof_neighbors);

  }

  static PyObject* flcbdfWrappersBuildDiscontinuousGalerkinLocal2GlobalMappings(PyObject* self,
                                                                                PyObject* args)
  {
    using namespace std;
    int nDOF_element;
    int nDOF_all_processes=0, nDOF_subdomain=0, max_dof_neighbors=0;
    PyObject *cmesh,*subdomain_cmesh,
      *elementOffsets_subdomain_owned,
      *elementNumbering_subdomain2global,
      *dg_dof_offsets_subdomain_owned,
      *dgNumbering_subdomain2global,
      *dg_subdomain_l2g;
    if (!PyArg_ParseTuple(args,
                          "iOOOOOOO",
                          &nDOF_element,
                          &cmesh,
                          &subdomain_cmesh,
                          &elementOffsets_subdomain_owned,
                          &elementNumbering_subdomain2global,
                          &dg_dof_offsets_subdomain_owned,
                          &dg_subdomain_l2g,
                          &dgNumbering_subdomain2global))

      return NULL;
    buildDiscontinuousGalerkinSubdomain2GlobalMappings(PROTEUS_COMM_WORLD,
                                                       MESH(cmesh),
                                                       IDATA(elementOffsets_subdomain_owned),
                                                       IDATA(elementNumbering_subdomain2global),
                                                       nDOF_element,
                                                       nDOF_all_processes,
                                                       nDOF_subdomain,
                                                       max_dof_neighbors,
                                                       IDATA(dg_dof_offsets_subdomain_owned),
                                                       IDATA(dg_subdomain_l2g),
                                                       IDATA(dgNumbering_subdomain2global));


    return Py_BuildValue("iii",
                         nDOF_all_processes,
                         nDOF_subdomain,
                         max_dof_neighbors);

  }

  static PyMethodDef flcbdfWrappersOldMethods[] = {
    { "globalSum",
      flcbdfWrappersGlobalSum,
      METH_VARARGS,
      "sum the value over all subdomains(processes)"},
    { "globalMax",
      flcbdfWrappersGlobalMax,
      METH_VARARGS,
      "take the max of the value over all subdomains(processes)"},
    { "globalMin",
      flcbdfWrappersGlobalMin,
      METH_VARARGS,
      "take the max of the value over all subdomains(processes)"},
    { "partitionElements",
      flcbdfWrappersPartitionElements,
      METH_VARARGS,
      "partition the mesh using an element-based partitioning"},
    { "partitionNodes",
      flcbdfWrappersPartitionNodes,
      METH_VARARGS,
      "partition the mesh using a node-based partitioning"},
    {"convertPUMIPartitionToPython",
     flcbdfWrappersConvertPUMIPartitionToPython,
     METH_VARARGS,
     "Convert C structures to python for PUMI partitioned mesh"},
    { "partitionNodesFromTetgenFiles",
      flcbdfWrappersPartitionNodesFromTetgenFiles,
      METH_VARARGS,
      "partition the mesh using a node-based partitioning"},
    { "buildQuadraticLocal2GlobalMappings",
      flcbdfWrappersBuildQuadraticLocal2GlobalMappings,
      METH_VARARGS,
      "create quadratic C0 finite element local to global mapping and subdomain 2 global mapping"},
    { "buildQuadraticCubeLocal2GlobalMappings",
      flcbdfWrappersBuildQuadraticCubeLocal2GlobalMappings,
      METH_VARARGS,
      "create quadratic C0 finite element local to global mapping and subdomain 2 global mapping for cubes"},
    { "buildDiscontinuousGalerkinLocal2GlobalMappings",
      flcbdfWrappersBuildDiscontinuousGalerkinLocal2GlobalMappings,
      METH_VARARGS,
      "create quadratic C0 finite element local to global mapping and subdomain 2 global mapping"},
    { NULL,NULL,0,NULL}
  };

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
  PyMODINIT_FUNC
  initflcbdfWrappersOld(void)
  {
    PyObject *m,*d,*c_api_object;
    static void* PyFLCBDFWrappersOld_API[1];
    m = Py_InitModule3("flcbdfWrappersOld",
                       flcbdfWrappersOldMethods,
                       "flcbdf wrappers module");
    d = PyModule_GetDict(m);
    import_array();
    // ensure PETSc, then DAETK, are initialized
    // PETSc first, via the proteus.Comm module
    //this wasn't working anyway...
    //PyRun_SimpleString("proteus.Comm.init()");

    // Set up default Proteus communicator
    PROTEUS_COMM_WORLD = PETSC_COMM_WORLD;

    PyFLCBDFWrappersOld_API[0] = (void*)(&PROTEUS_COMM_WORLD);
    c_api_object = PyCObject_FromVoidPtr((void*)PyFLCBDFWrappersOld_API,NULL);
    PyModule_AddObject(m,"_C_API",c_api_object);
  }
}
/** @} */
