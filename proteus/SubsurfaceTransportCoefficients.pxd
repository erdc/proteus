# A type of -*- python -*- file
cdef extern from "SubsurfaceTransportCoefficients.h":
    cdef void piecewiseLinearTableLookup(double x,
				         int nv,
				         int* start,
				         double* y,
				         double* dy,
				         const double* xv,
				         const double* yv)

    cdef int calculateRusanovFluxSaturationEquationIncomp_PWC(double safetyFactor,
						              int nSpace,
						              #physical information
						              #psk type
						              int pskModelFlag,
						              int nParams,
						              const int* rowptr,
						              const int* colind,
						              const int* materialTypes,
						              double muw,
						              double mun,
						              const double* omega,
						              const double* Kbar,
						              double b,
						              const double* rwork_psk,
						              const double* rwork_psk_tol,
						              const double* rwork_density_w,
						              const double* rwork_density_n,
						              const double* g,
						              #velocity at interface 
						              const double* ebq_global_qt,
						              #bounds on characteristic speed for elements
						              const double* q_lambda_bar,
						              #mesh information
						              int nElements_global,
						              int nElementBoundaries_element,
						              int nInteriorElementBoundaries_global,
						              int nExteriorElementBoundaries_global,
						              int nQuadraturePoints_element, 
						              int nQuadraturePoints_elementBoundary,
						              const int* interiorElementBoundaries,
						              const int* exteriorElementBoundaries,
						              const int* elementBoundaryElements,
						              const int* elementBoundaryLocalElementBoundaries,
						              const double* n,
						              #solution information
						              const double * q_u,
						                #element quadrature
						              const double* q_dV,
						                #int nDOF_trial_element,
						              #const int* u_l2g,
						              #const double* u_dof,
						              #boundary information
						              const int * isDOFBoundary,
						              const double* bc_u,
						                #output
						              double* flux)
