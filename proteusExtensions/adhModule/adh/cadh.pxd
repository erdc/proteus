cdef extern from "cadhimpl.h":
    ctypedef struct ADH_Input:
        #only declare date members that we need access to
        pass
    ctypedef struct ADH_PETSc_Interface:
        #
        pass
    ctypedef struct ADH_NumericalSolution:
        #only declare date members that we need access to
        pass
    int cadhMain(int argc, char* argv[])
    void readBC(ADH_Input*, char* filebase)
    void messg_finalize_proteus()
    int ADH_NumericalSolution_step(ADH_NumericalSolution* self)
    int ADH_NumericalSolution_stepTaken(ADH_NumericalSolution* self)
    int ADH_NumericalSolution_init(ADH_NumericalSolution* self, int argc, char* argv[])
    int ADH_NumericalSolution_dealloc(ADH_NumericalSolution* self)
    int ADH_NumericalSolution_calculateSolution(ADH_NumericalSolution* self)

    ctypedef struct ADH_OneLevelTransport:
        #Newton increment
        double *sol
        #Nonlinear residual
        double *residual
        #base size of system (not including multiple dofs per node)
        int isize
    int realloc_fe_main(ADH_OneLevelTransport* fe)
    int ADH_OneLevelTransport_init(ADH_OneLevelTransport* self)
    int ADH_OneLevelTransport_dealloc(ADH_OneLevelTransport* self)
    #Call the ADH Newton solves for all the models. This will kick out
    #with a failure for any model so the caller needs handle nonlinear
    #solver failures
    int ADH_OneLevelTransport_solve(ADH_OneLevelTransport* self)

    #solution details
    int get_ADH_NumberOfComponents()
    int get_ADH_NumberOfDegreesOfFreedom()
    double get_ADH_newTimeLevel()
    double get_ADH_prevTimeLevel()
    double get_ADH_dt()
    double get_ADH_tinitial()
    #pointer to a solution component, depends on the model
    void* get_ADH_solutionComponent(int ci)
    #logical dimensionality of solution component ci 
    int get_ADH_solutionComponentDimension(int ci)
    #starting component number for vector-valued quantity that contains component ci
    int get_ADH_firstSolutionComponentForVectorComponent(int ci)
    int get_ADH_ADAPTION_flag()
    
    #mesh info
    int get_ADH_nElements_global()
    int get_ADH_nNodes_global()
    int get_ADH_nSpace_global()
    int get_ADH_nNodes_element()
    int get_nodeArray(double* nodeArray)
    int get_elementNodesArray(int* elementNodesArray)
    int get_elementMaterialTypes(int* elementMaterialTypes)
    int get_nodeMaterialTypes(int* nodeMaterialTypes)

    #ADH_PETSc_Interface
    int ADH_PETSc_Interface_init(ADH_PETSc_Interface* self, ADH_OneLevelTransport* cadh_transport)
#     void ADH_PETSC_matrix_resize(ADH_PETSc_Interface* self,
#                                  ADH_OneLevelTransport* cadh_transport,
#                                  KSP* ksp, 
#                                  Mat* matrix,
#                                  int nsys, 
#                                  Vec* sol, 
#                                  Vec* residual);
#     void ADH_PETSc_matrix_destroy(KSP* ksp, Mat *matrix, Vec* sol, Vec* residual)
    int ADH_PETSc_Interface_update(ADH_PETSc_Interface* self, ADH_OneLevelTransport* cadh_transport)
    int ADH_PETSc_updateStorage(ADH_PETSc_Interface* self, ADH_OneLevelTransport* cadh_transport)
    int ADH_PETSc_Interface_dealloc(ADH_PETSc_Interface* self)
    int ADH_OneLevelTransport_updateStorage(ADH_OneLevelTransport* self)
    
    
