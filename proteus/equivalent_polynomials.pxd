# A type of -*- python -*- file
from libcpp cimport bool
cdef extern from "equivalent_polynomials.h" namespace "equivalent_polynomials":
    cdef cppclass cSimplex "equivalent_polynomials::Simplex"[nSpace,nP,nQ]:
      cSimplex "Simplex"()
      void calculate(double* phi_dof, double* phi_nodes, double* xi_r);
      double* get_H()
      double* get_ImH()
      double* get_D()
      bool inside_out
