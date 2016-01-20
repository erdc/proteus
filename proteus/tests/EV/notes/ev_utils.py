#! /usr/bin/env python
"""
store some routines for trying out an EV approximation in 1d
"""

import numpy as np
import scipy
from scipy import linalg,sparse
from scipy.sparse import linalg


def choose_stable_dt(uelem,dflux,he,x,t):
    """
    pick time steps according to a target CFL
    """
    wave_speeds = dflux(uelem,x,t)
    wave_speeds = np.max(np.absolute(wave_speeds),1)/he
    fastest = wave_speeds.max()
    return 1.0/(fastest+ 1.0e-10)



"""
### SSP RK3 update
\begin{align}
\mathbf{k}_1 &= \mathbf{F}_h(\mathbf{U}^n,t^n) \\
\mathbf{k}_2 &= \mathbf{F}_h(\mathbf{U}^n + a_{21}\Delta t\mathbf{k}_1,t^n+ c_2\Delta t) \\
\mathbf{k}_3 &= \mathbf{F}_h(\mathbf{U}^n + a_{31}\Delta t\mathbf{k}_1 + a_{32}\Delta t\mathbf{k}_2,t^n+ c_3\Delta t) \\
\mathbf{U}^{n+1} &= \mathbf{U}^n+\Delta t\sum_{i=1}^3b_i\mathbf{k}_i
\end{align}
with
\begin{eqnarray}
c_2 =1,  &a_{21} = 1 &\\
c_3 =1/2,& a_{31} = 1/4, & a_{32} = 1/4 \\
b_1 =1/6,& b_2 = 1/6, & b_3= 2/3
\end{eqnarray}
"""

def ssp_rk1(un,tn,dt,rhs,constraint_ind,constraint_vals):
    k1 = rhs(un,tn)
    u1 = un+dt*k1; 
    u1[constraint_ind]=constraint_vals
    
    return u1,tn+dt



def ssp_rk3(un,tn,dt,rhs,constraint_ind,constraint_vals):
    c2=1.; a21=1.;
    c3=0.5; a31=0.25; a32=0.25
    b1=1.0/6; b2=1.0/6.0; b3=2.0/3.0
    
    k1 = rhs(un,tn)
    u1 = un+a21*dt*k1; 
    u1[constraint_ind]=constraint_vals#should update these in time
    
    k2 = rhs(u1,tn+c2*dt)
    u2 = un + a31*dt*k1 + a32*dt*k2
    u2[constraint_ind]=constraint_vals
    
    k3 = rhs(u2,tn + c3*dt)
    
    uh = un + dt*(b1*k1+b2*k2+b3*k3)
    uh[constraint_ind]=constraint_vals
    
    return uh,tn+dt


def build_graph_laplacian(he,Ndof,nu):
    """
    \begin{align}
                            & -\frac{1}{n_e}|\Omega_e|, \ i \ne j, \ i,j \in \mathcal{I}(\Omega_e) \\
    b_{e}(w_{h,j},w_{h,i}) =& \ |\Omega_e|, \  i = j, \ i,j \in \mathcal{I}(\Omega_e) \\
                            & \ 0, \ \mbox{ otherwise}
    \end{align}

    """
    graph_diag=np.zeros((Ndof,),'d')
    graph_diag[1:-1]=nu[0:-1]*he[0:-1]+nu[1:]*he[1:] #diagonal interior
    graph_diag[0]=nu[0]*he[0]; graph_diag[-1]=nu[-1]*he[-1]#account for boundary
    diagonals=[graph_diag,-nu*he,-nu*he]; offsets=[0,-1,1]
    B = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    return B

def low_order_divergence(un,x,tn,he,l2g,K,C,Minv,out_ind,out_normals,
                         flux,dflux,low_visc,stab_operator,
                         include_inverse_mass=True):
    """
    compute low-order spatial divergence terms (integrated by parts)
    This is FI, I=0,\dots,Nn-1 in Guermond
    If include_inverse_mass is true, multiply by (M^L)-1
    Also returns low-order viscosity
    """
    Ndof = un.shape[0]
    #explicit diffusion
    b = K.dot(un)

    ##explicit advection
    fn = flux(un,x,tn)
    b -= C.dot(fn)
    ##outflow
    b[out_ind] += fn[out_ind]*out_normals
    ##low-order artificial viscosity
    dfn=dflux(un,x,tn)
    nu_L = low_visc(he,l2g,dfn)
    #low-order stabilization
    B_L = stab_operator(he,Ndof,nu_L)
    b  += B_L.dot(un)
    if include_inverse_mass:
        b   = Minv.dot(b)
    return b,nu_L

def low_order_step(un,x,tn,dt,he,l2g,K,C,Minv,dir_ind,dir_values,out_ind,out_normals,
                   flux,dflux,low_visc,stab_operator):
    ##divergence
    F,nu_L = low_order_divergence(un,x,tn,he,l2g,K,C,Minv,out_ind,out_normals,
                                  flux,dflux,low_visc,stab_operator,include_inverse_mass=True)
    #add previous solution
    uhL = un - dt*F
    ##enforce constraints
    uhL[dir_ind]=dir_values

    return uhL,nu_L


def high_order_divergence(un,nu_L,x,tn,R_e,eta,he,l2g,K,C,MinvH,out_ind,out_normals,
                          flux,dflux,low_visc,high_visc,edge_jump,entropy_normalization,
                          stab_operator,
                          include_inverse_mass=True):
    """
    compute high-order spatial divergence terms (integrated by parts)
    This is G_I, I=0,\dots,Nn-1 in Guermond
    If include_inverse_mass is true, multiply by (M^H)-1
    Also returns high-order viscosity
    """
    Ndof = un.shape[0]

    #explicit diffusion
    b = K.dot(un)

    ##explicit advection
    fn = flux(un,x,tn)
    b -= C.dot(fn)
    ##outflow
    b[out_ind] += fn[out_ind]*out_normals

    ##high-order artificial viscosity
    dfn=dflux(un,x,tn)
    Re_norm = np.amax(np.absolute(R_e),1)
    eta_jump = edge_jump(he,l2g,eta,dfn)
    eta_diff = entropy_normalization(eta)

    Jf = np.amax(np.absolute(eta_jump[l2g]),1)
    nu_H = high_visc(nu_L,Re_norm,eta_diff,Jf)

    #high-order stabilization
    B_H = stab_operator(he,Ndof,nu_H)
    b  += B_H.dot(un)
 
    if include_inverse_mass:
        b = MinvH.dot(b)
    return b,nu_H

def high_order_step(un,nu_L,x,tn,dt,R_e,eta,he,l2g,K,C,Minv,dir_ind,dir_values,out_ind,out_normals,
                    flux,dflux,low_visc,high_visc,stab_operator):
    G,nu_H = high_order_divergence(un,nu_L,x,tn,R_e,eta,he,l2g,K,C,MinvH,out_ind,out_normals,
                                   flux,dflux,low_visc,high_visc,stab_operator,
                                   include_inverse_mass=True)
    #add previous solution
    uhH  = un.copy(); uhH -= dt*G
    
    ##enforce constraints
    uhH[dir_ind]=dir_values

    return uhH,nu_H


def build_linear_operators(Ndof,a0,he,l2g,dir_ind,out_ind,out_normals):
    main_diag=np.ones((Ndof,)); offdiag=np.ones((Ndof-1),)
    lower_diag=offdiag; upper_diag=offdiag
    ##remove strong boundary conditions from entries
    main_diag_bc=main_diag.copy()
    main_diag_bc[dir_ind]=0.0
    upper_diag_bc=offdiag.copy()
    lower_diag_bc=offdiag.copy()
    for dir in dir_ind:
        if dir < Ndof-1:
            upper_diag_bc[dir]=0.
        if dir > 0:
            lower_diag_bc[dir-1]=0.
    ##account for outflow boundary
    main_diag_bc[out_ind] *= 0.5
    adv_main_diag=np.zeros(main_diag.shape,'d')
    adv_main_diag[out_ind]=1.0
    ##Stiffness
    diagonals=[2.0*a0/he[0]*main_diag_bc,-a0/he[0]*lower_diag_bc,-a0/he[0]*upper_diag_bc]; offsets=[0,-1,1]
    K = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##mass
    diagonals=[2.0/3.0*he[0]*main_diag,1.0/6.0*he[0]*lower_diag,1.0/6.0*he[0]*upper_diag]; offsets=[0,-1,1]
    Mc = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##lumped mass matrix 
    diagonals=[he[0]*main_diag]; offsets=[0]
    diagonals[0][0]=5.0/6.0*he[0]; diagonals[0][-1]=5.0/6.0*he[0]
    M = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    diagonals[0]=1.0/diagonals[0]
    Minv = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##Approximate inverse for the consistent mass matrix
    MB=M-Mc
    MB = MB.dot(Minv)
    MBplusEye = MB + scipy.sparse.eye(Ndof)
    MinvH = Minv.dot(MBplusEye)

    ##advection
    diagonals=[0.5*adv_main_diag,0.5*lower_diag_bc,-0.5*upper_diag_bc]; offsets=[0,-1,1]
    C = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')

    ##Dirichlet
    dir_dummy=np.ones(dir_ind.shape,'d')
    D = scipy.sparse.coo_matrix((dir_dummy,(dir_ind,dir_ind)),shape=(Ndof,Ndof)).tocsr()

    return Mc,M,Minv,MinvH,MB,C,K,D

def build_linear_operators_no_dir(Ndof,a0,he,l2g,dir_ind,out_ind,out_normals):
    
    main_diag=np.ones((Ndof,)); offdiag=np.ones((Ndof-1),)
    upper_diag=offdiag; lower_diag=offdiag;
    ##remove strong boundary conditions from entries (may not need anymore)
    main_diag_bc=main_diag.copy()
    main_diag_bc[dir_ind]=0.0
    upper_diag_bc=offdiag.copy()
    lower_diag_bc=offdiag.copy()
    for dir in dir_ind:
        if dir < Ndof-1:
            upper_diag_bc[dir]=0.
        if dir > 0:
            lower_diag_bc[dir-1]=0.
            ##account for outflow boundary
    main_diag_bc[out_ind] *= 0.5
    adv_main_diag=np.zeros(main_diag.shape,'d')
    adv_main_diag[out_ind]= out_normals
    #not removing dirichlet boundaries from matrices
    adv_main_diag[dir_ind]=-1.0
    ##Stiffness
    diagonals=[2.0*a0/he[0]*main_diag,-a0/he[0]*lower_diag,-a0/he[0]*upper_diag]; offsets=[0,-1,1]
    K = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##consistent mass
    diagonals=[2.0/3.0*he[0]*main_diag,1.0/6.0*he[0]*lower_diag,1.0/6.0*he[0]*upper_diag]; offsets=[0,-1,1]
    Mc = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##lumped mass matrix 
    diagonals=[he[0]*main_diag]; offsets=[0]
    diagonals[0][0]=5.0/6.0*he[0]; diagonals[0][-1]=5.0/6.0*he[0]
    M = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    diagonals[0]=1.0/diagonals[0]
    Minv = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##Approximate inverse for the consistent mass matrix
    MB=M-Mc
    MB = MB.dot(Minv)
    MBplusEye = MB + scipy.sparse.eye(Ndof)
    MinvH = Minv.dot(MBplusEye)

    ##advection
    diagonals=[0.5*adv_main_diag,0.5*lower_diag,-0.5*upper_diag]; offsets=[0,-1,1]
    C = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    ##Dirichlet
    dir_dummy=np.ones(dir_ind.shape,'d')
    D = scipy.sparse.coo_matrix((dir_dummy,(dir_ind,dir_ind)),shape=(Ndof,Ndof)).tocsr()

    return Mc,M,Minv,MinvH,MB,C,K,D




def build_A1_flux_matrix(U,nu_H,nu_L,he,dt):
    Ndof = U.shape[0]
    #A1
    diagonals=[-dt*(nu_L-nu_H)*(U[1:]-U[0:-1])*he,-dt*(nu_L-nu_H)*(U[0:-1]-U[1:])*he]
    offsets=[1,-1]
    A1 = scipy.sparse.diags(diagonals,offsets,shape=(Ndof,Ndof),format='csr')
    return A1

def build_A2_flux_matrix_new(dt,B,G,A2):
    Ndof = G.shape[0]
    for I in range(Ndof):
        for J in range(max(0,I-1),min(Ndof,I+2)):
            #print 'setting I={0} J={1}'.format(I,J)
            A2[I,J]=-dt*(B[I,J]*G[J]-B[J,I]*G[I])


    return A2
def build_A1_flux_matrix_new(U,nu_H,nu_L,he,dt,l2g,Ig,Jg):
    nDOF=U.shape[0]; 
    nN_local  = l2g.shape[1]
    #dense storage
    A1g = np.zeros((l2g.shape[0],nN_local**2),'d')
    for ii in range(nN_local):
        for jj in range(nN_local):
            A1g[:,jj*nN_local+ii] = -dt*(nu_L-nu_H)*(U[l2g][:,jj]-U[l2g][:,ii])*he
    return scipy.sparse.coo_matrix((A1g.flat,(Ig.flat,Jg.flat)),shape=(nDOF,nDOF)).tocsr()

"""
Write SSP RK schemes as a series of Forward Euler steps

RK3:
\begin{align}
u^0 &=u^n, t^0=  t^n \\
u^1 &=u^0 + \Delta t L(u^0,t^0) \rightarrow  \mbox{ \texttt{fe\_step}}(u^0,t^0,\Delta t) \\
\tilde{u}^2 &= u^1 + \Delta t L(u^1,t+\Delta t)  \rightarrow  \mbox{ \texttt{fe\_step}}(u^1,t^0+\Delta t,\Delta t) \\
u^2 &= \frac{3}{4}u^0 + \frac{1}{4}\tilde{u}^2 \\
\tilde{u}^3 &= u^2 + \Delta t L(u^2,t+\frac{\Delta t}{2})  \rightarrow  \mbox{ \texttt{fe\_step}}(u^2,t^0+\Delta t/2,\Delta t) \\
u^3 &= \frac{1}{3}u^0 + \frac{2}{3}\tilde{u}^3 \\
u^{n+1} &= u^3
\end{align}
"""

def ssp_rk1_ev(un,tn,dt,fe_step,R_e,eta,entropy_residual,he,l2g,x,entropy,dflux):
    u1 = fe_step(un,tn,dt,R_e,eta)

    R_e = entropy_residual(he,l2g,dt,tn,x,u1,un,entropy,dflux)
    eta = entropy(u1,x,tn)


    return u1,tn+dt,R_e,eta

def ssp_rk3_ev(un,tn,dt,fe_step,R_e,eta,entropy_residual,he,l2g,x,entropy,dflux):
    """
    RK3:
    \begin{align}
    u^0 &=u^n, t^0=  t^n \\
    u^1 &=u^0 + \Delta t L(u^0,t^0) \rightarrow  \mbox{ \texttt{fe\_step}}(u^0,t^0,\Delta t) \\
    \tilde{u}^2 &= u^1 + \Delta t L(u^1,t+\Delta t)  \rightarrow  \mbox{ \texttt{fe\_step}}(u^1,t^0+\Delta t,\Delta t) \\
    u^2 &= \frac{3}{4}u^0 + \frac{1}{4}\tilde{u}^2 \\
    \tilde{u}^3 &= u^2 + \Delta t L(u^2,t+\frac{\Delta t}{2})  \rightarrow  \mbox{ \texttt{fe\_step}}(u^2,t^0+\Delta t/2,\Delta t) \\
    u^3 &= \frac{1}{3}u^0 + \frac{2}{3}\tilde{u}^3 \\
    u^{n+1} &= u^3
    \end{align}
    """
    one_third=1.0/3.0 ; two_third=2.0/3.0
    u1 = fe_step(un,tn,dt,R_e,eta)

    R_e = entropy_residual(he,l2g,dt,tn,x,u1,un,entropy,dflux)
    eta = entropy(u1,x,tn)

    ut2= fe_step(u1,tn+dt,dt,R_e,eta)
    u2 = 0.75*un+0.25*ut2

    R_e = entropy_residual(he,l2g,dt,tn,x,u2,u1,entropy,dflux)
    eta = entropy(u2,x,tn)
    
    ut3= fe_step(u2,tn+0.5*dt,dt,R_e,eta)
    u3 = one_third*un + two_third*ut3

    R_e = entropy_residual(he,l2g,dt,tn,x,u3,u2,entropy,dflux)
    eta = entropy(u2,x,tn)

    return u3,tn+dt,R_e,eta
