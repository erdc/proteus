cimport numpy as np
cdef extern from "math.h":
    double fmax(double a, double b)
    double fmin(double a, double b)

def build_A2_flux_matrix(double dt,np.ndarray G,np.ndarray row,np.ndarray col,np.ndarray BIJ,np.ndarray AIJ):
    """
    Assemble flux difference for matrix. Assumes Identical layout for B and A
    """
    cdef int I,m,J,mm,II
    cdef int Nn=row.shape[0]-1
    cdef int* rowptr=<int*> row.data
    cdef int* colptr=<int*> col.data
    cdef double* Gptr=<double*> G.data
    cdef double* BIJptr=<double*> BIJ.data
    cdef double* AIJptr=<double*> AIJ.data

    for I in range(Nn):
        for m in range(rowptr[I],rowptr[I+1]):
            J = colptr[m]
            AIJptr[m] = -dt*BIJptr[m]*G[J]
            #mwf debug
            #print '({0},{1}) BIJ={2},G[{1}]={3}'.format(I,J,BIJptr[m],G[J])
            #find BJI
            for mm in range(rowptr[J],rowptr[J+1]):
                II = colptr[mm]
                if II==I:
                    AIJptr[m] += dt*BIJptr[mm]*G[I]
                    #mwf debug
                    #print '({0},{1}) --> ({2},{3} BJI={4},G[{5}]={6}'.format(I,J,J,II,BIJptr[mm],I,G[I])
                    break
    #
    return AIJ

def apply_fct(np.ndarray Un, np.ndarray UL, np.ndarray row, np.ndarray col, np.ndarray MIJ, np.ndarray  AIJ, 
              np.ndarray Rm, np.ndarray Rp):
    """
    Attempt FCT algorithm for Anti-diffusive fluxes stored in matrix A assuming CSR format.
    Assumes M is diagonal
    """
    cdef int* rowptr=<int*> row.data #pointer to nonzero entries for row I
    cdef int* colptr=<int*> col.data #column indices
    cdef double* MIJptr=<double*> MIJ.data
    cdef double* AIJptr=<double*> AIJ.data
    cdef double* Rmptr =<double*> Rm.data
    cdef double* Rpptr =<double*> Rp.data
    
    #local variables 
    cdef int I,m,J
    cdef int Nn=row.shape[0]-1
    cdef double eps = 1.0e-12
    
    cdef double U_I,U_J,UL_I,U_max,U_min,QIp,QIm,PIp,PIm,MI,RIp,RIm,LIJ
    for I in range(Nn):
        #1. Calculate local max and min
        U_I=Un[I]; U_max=U_I; U_min=U_I
        for m in range(rowptr[I],rowptr[I+1]):
            J = colptr[m]
            U_max = fmax(U_max,Un[J])
            U_min = fmin(U_min,Un[J])
        #2. Compute distance to local extrema
        MI = MIJptr[I]; UL_I=UL[I]
        #mwf check this
        QIp = MI*(U_max-UL_I); QIm = MI*(U_min - UL_I)
        
        #3. Compute sums of positive and negative anti-diffusive fluxes
        PIp=0.; PIm=0.;
        for m in range(rowptr[I],rowptr[I+1]):
            PIp += fmax(0.,AIJptr[m])
            PIm += fmin(0.,AIJptr[m])
        #4. Compute the nodal correction factors
        RIp=1.0; RIm=1.0
        if PIp > eps: 
            RIp = fmin(1.0,QIp/PIp)
        if PIm < -eps:
            RIm = fmin(1.0,QIm/PIm)
        #save for limiting step
        Rpptr[I]=RIp; Rmptr[I]=RIm
        #if QIp > eps:
        #    print "I={0}, U_max={1}, UL_I={2}, QIp={3}, PIp={4}, RIp={5}".format(I,U_max,UL_I,QIp,
        #                                                                   PIp,RIp)
        #if QIm < -eps:
        #    print "I={0}, U_min={1}, UL_I={2}, QIm={3}, PIm={4}, RIm={5}".format(I,U_min,UL_I,QIm,
        #                                                                         PIm,RIm)
        #assert QIp >= 0.0, "I={0}, U_max={1}, UL_I={2}".format(I,U_max,UL_I)
        #assert QIm <= 0.0, "I={0}, U_min={1}, UL_I={2}".format(I,U_min,UL_I)

        
    #end first loop
    
    ##compute limiting factor
    for I in range(Nn):
        for m in range(rowptr[I],rowptr[I+1]):
            J = colptr[m]
            LIJ = 1.0
            if AIJptr[m] >= 0.0:
                LIJ = fmin(Rpptr[I],Rmptr[J])
            else:
                LIJ = fmin(Rmptr[I],Rpptr[J])
            AIJptr[m] *= LIJ
    #
    return Rp,Rm

def build_div_representation(np.ndarray row, np.ndarray col):
    """
    Given sparsity information for relating scalar spaces in operator like 
    (phi_i,phi_j)_\Omega, phi_i \in W, phi_j \in V
    compute sparsity information for a 'divergence' matrix
    (phi_i,grad phi_j)_Omega 

    """
    pass

def build_graph_laplacian_element(np.ndarray q_dV, np.ndarray b_e, int nDOF_trial_element):
    """
    Build element contribution for Graph laplacian
    \begin{align}
                            & -\frac{1}{n_e}|\Omega_e|, \ i \ne j, \ i,j \in \mathcal{I}(\Omega_e) \\
    b_{e}(w_{h,j},w_{h,i}) =& \ |\Omega_e|, \  i = j, \ i,j \in \mathcal{I}(\Omega_e) \\
                            & \ 0, \ \mbox{ otherwise}
    \end{align}

    """
    cdef int i,j,iq,eN,nElements_global=q_dV.shape[0], nQuadrature_element=q_dV.shape[1]
    cdef double vol,bij,bjj,ne_inv
    assert b_e.shape[0] == nElements_global
    assert b_e.shape[1] == nDOF_trial_element
    assert b_e.shape[2] == nDOF_trial_element
    
    cdef double* dVptr=<double*> q_dV.data
    cdef double* b_eptr=<double*> b_e.data
    ne_inv = 1.0/(nDOF_trial_element - 1.0)
    
    for eN in range(nElements_global):
        vol = 0.0
        for iq in range(nQuadrature_element):
            vol += dVptr[eN*nQuadrature_element+iq]
        bji = -vol*ne_inv
        bjj = vol
        for j in range(nDOF_trial_element):
            for i in range(nDOF_trial_element):
                b_eptr[eN*nDOF_trial_element*nDOF_trial_element + j*nDOF_trial_element + i] = bji
            #overwrite diagonal
            b_eptr[eN*nDOF_trial_element*nDOF_trial_element + j*nDOF_trial_element + j] = bjj
    return b_e

