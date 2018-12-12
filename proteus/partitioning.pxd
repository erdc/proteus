cimport mesh as cppm

cdef extern from "partitioning.h":
    cdef int partitionNodes(cppm.Mesh& mesh,
                            int nNodes_overlap)
