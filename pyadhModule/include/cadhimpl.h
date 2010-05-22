#ifndef CADH_INCLUDE_H
#define CADH_INCLUDE_H

#ifndef ADH_SHARE_INCLUDED
#include "share_extern.h"
#endif
#include "petscksp.h"

/** The C representation of of ADH_Input, a class for
    representing just the input, which may consist
    of multiple equation sets. This is basically the local variables from adh.c */
typedef struct
{
  /*char binfile[MAXLINE]; */
  int iout;			/* output counter */
  int ixy;			/* counter for the XY series */
  int irn;			/* run name counter */
  int interval;		/* interval var */
  int hot_count;		/* keeps track of when a hot start file needs to be printed */
  FILE *fp_in;		/* input file */
  FILE *fp_sflx_in;	/* input file */
  FILE *fp_sflx_out;	/* output surface heat flux */
#ifdef _ADH_XT3
  char *stdbuf;
#endif
} ADH_Input;

typedef struct
{
  int *nnodeLocal;
  int *nnodeLocal_old;
  int *globalNodeNumberStart;
  int *globalNodeNumber;
  int mappingHasChanged;
  ISLocalToGlobalMapping localToGlobalMapping;
  Mat L[10];
  Vec du[10],r[10];
  KSP ksp[10];
} ADH_PETSc_Interface;

/** The C representation of of ADH_NumericalSolution, a class for
    representing the entire ADH numerical solution, which may consist
    of multiple equation sets. This is basically the local variables from adh.c */
typedef struct
{
  /*char binfile[MAXLINE]; */
  int iout;			/* output counter */
  int ixy;			/* counter for the XY series */
  int irn;			/* run name counter */
  int interval;		/* interval var */
  int hot_count;		/* keeps track of when a hot start file needs to be printed */
  FILE *fp_in;		/* input file */
  FILE *fp_sflx_in;	/* input file */
  FILE *fp_sflx_out;	/* output surface heat flux */
#ifdef _ADH_XT3
  char *stdbuf;
#endif
} ADH_NumericalSolution;


extern int cadhMain(int argc, char* argv[]);
extern void readBC(ADH_Input* self, char* filebase);
extern void messg_finalize_pyadh(void);
extern int ADH_NumericalSolution_step(ADH_NumericalSolution* self);
extern int ADH_NumericalSolution_stepTaken(ADH_NumericalSolution* self);
extern int ADH_NumericalSolution_init(ADH_NumericalSolution* self, int argc, char* argv[]);
extern int ADH_NumericalSolution_dealloc(ADH_NumericalSolution* self);
extern int ADH_NumericalSolution_calculateSolution(ADH_NumericalSolution* self);



/**The C representation of an ADH_OneLevelTransport object. These are
   basically the variables in fe_main, except that I added a pointer
   to a cADH_NumericalSolution object to make sure none of the global
   variables get deleted before where done with the transport
   object. */ 
typedef struct
{
  int isize;/* the size of the arrays */
  SPARSE_VECT *matrix;/* the matrix */
  double *diagonal;/* the diagonal */
  int *bc_mask;/* this indicates a dirichlet boundary condition */
  double *nodal_production;/* the turbulent prodution for NS distributed to the nodes */
  double *nodal_viscosity;/* the nodal viscosity */
  double *nodal_volume;/* the volume of the nodes */
  double *residual;/* the residual - which is also the right hand side */
  double *sol;/* the solution */
  double *scale_vect;/* the scale vector for the jacobi preconditioning */
  VECT **bed_load_3d;/* the total bed load for 3 dimensions */
  double *ustar;
  double *ald;/* temporary array to store active layer distribution in the form that the sed */
  double *asd;/* temporary array to store active stratum distribution in the form that the sed */
  double *blt;/* temporary array to store the bed layer thickness in the form that the sed */
  double *bl_por;
  double *bl_bd;
  double *bl_sg;
  double *bl_ces;
  double *bl_erc;
  double *bl_ere;
  double **bld;/* temporary array to store the bed layer distribution in the form that the sed */
  double **cbp;/* temporary array to store the consolidation bed properties in the form that the sed */
  double *sum_snd_ald;/* used to determine the sum of sand fractions in the active layer */
  int *bed_character;/* used to determine if the active window will be considered clay or sand */
  int *aldkf;
  int *asdkf;
  double *sum_nodal_source;/* used to determine the total source for sands combined */
  double *sum_snd_con; /* used to determine the total source for sands combined */
  double **source_coef;/* the equilibrium concentration for sands, a source for clays */
  double **prob_dep;/* the probability of deposition for clays */
  VECT2D **eq_bed_load;
  double **bedload_eq_con;
  double **bl_thick;
  double **bl_slc;
  VECT2D **bl_vel;
  double *sum_blt;
  int *bed_flag;/*flag set to 1 if the node lies on the bed */
  double *old_sed_displacement;/* the displacement from the previous sed timestep */
  double *old_sed_active_stratum_ceiling;/* the ASC from the previous sed timestep */
  double **old_sed_active_layer_distribution;/* the ALD from the previous sed timestep */
  double **old_sed_active_stratum_distribution;/* the ASD from the previous sed timestep */
  double **old_it_ald;	/* the active layer distribution from th eprevious iteration */
  double **old_it_rhs;	/* the rhs from th eprevious iteration */
  double **bl_eq_coef;
  double **critical_shear_vel;
  VECT2D *nodal_grad;
  double *elev_disp;
  double *fraction;
  double *old_fraction; 
  double *AL_thickness; 
  double *old_sed_ald;
  double *old_sed_active_layer_porosity;
  double *old_sed_active_layer_critical_erosion_shear;
  double *old_sed_active_layer_erosion_rate_constant;
  double *old_sed_active_layer_erosion_rate_exponent;
  double *sum_sed_decay;
  int *limiter;		/*the switch to determine if too much bed change is occurring */
  int *limiter2;		/*the switch to determine if bedsource should be included */
  double **old_it_concentration;	/* the concentration from the previous iteration */
  int *conc_converge;	/* counter to determine if concentration has converged */
} ADH_OneLevelTransport;

extern int realloc_fe_main(ADH_OneLevelTransport* fe);
extern int ADH_OneLevelTransport_init(ADH_OneLevelTransport* self);
extern int ADH_OneLevelTransport_dealloc(ADH_OneLevelTransport* self);
/** Call the ADH Newton solves for all the models. This will kick out
    with a failure for any model so the caller needs handle nonlinear
    solver failures */
extern int ADH_OneLevelTransport_solve(ADH_OneLevelTransport* self);

/*solution info*/
extern int get_ADH_NumberOfComponents(void);
extern int get_ADH_NumberOfDegreesOfFreedom(void);
extern double get_ADH_newTimeLevel(void);
extern double get_ADH_prevTimeLevel(void);
extern double get_ADH_dt(void);
extern double get_ADH_tinitial(void);
/*pointer to a solution component, depends on the model*/
extern void* get_ADH_solutionComponent(int ci);
/*logical dimensionality of solution component ci*/ 
extern int get_ADH_solutionComponentDimension(int ci);
/*starting component number for vector-valued quantity that contains component ci*/
extern int get_ADH_firstSolutionComponentForVectorComponent(int ci);
extern int get_ADH_ADAPTION_flag(void);
/*mesh info*/
extern int get_ADH_nElements_global(void);
extern int get_ADH_nNodes_global(void);
extern int get_ADH_nSpace_global(void);
extern int get_ADH_nNodes_element(void);
extern int get_nodeArray(double* nodeArray);
extern int get_elementNodesArray(int* elementNodesArray);
extern int get_elementMaterialTypes(int* elementMaterialTypes);
extern int get_nodeMaterialTypes(int* nodeMaterialTypes);
extern int ADH_PETSc_Interface_init(ADH_PETSc_Interface* self,
				   ADH_OneLevelTransport* cadh_transport);
extern void ADH_PETSc_matrix_resize(ADH_PETSc_Interface* self,
				    ADH_OneLevelTransport* cadh_transport,
				    KSP* ksp, 
				    Mat* matrix,
				    int nsys, 
				    Vec* sol, 
				    Vec* residual);
extern void ADH_PETSc_matrix_destroy(KSP* ksp, Mat *matrix, Vec* sol, Vec* residual);
extern int ADH_PETSc_updateStorage(ADH_PETSc_Interface* self,
				   ADH_OneLevelTransport* cadh_transport);
extern int ADH_PETSc_destroyStorage(ADH_PETSc_Interface* self);
extern int ADH_PETSc_Interface_update(ADH_PETSc_Interface* self,
				      ADH_OneLevelTransport* cadh_transport);
extern int ADH_PETSc_Interface_dealloc(ADH_PETSc_Interface* self);

#endif
