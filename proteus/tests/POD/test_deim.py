#!/usr/bin/env python

"""
some tests for initial DEIM implementation in proteus
"""
import numpy as np
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

from proteus import deim_utils

def get_burgers_ns(name,T=0.1,nDTout=10,archive_pod_res=None):
    import burgers_init as bu
    bu.physics.name=name
    bu.so.name = bu.physics.name
    #adjust default end time and number of output steps
    bu.T=T
    bu.nDTout=nDTout
    bu.DT=bu.T/float(bu.nDTout)
    bu.so.tnList = [i*bu.DT for i in range(bu.nDTout+1)]
    #request archiving of spatial residuals ...
    simFlagsList=None
    if archive_pod_res is not None:
        simFlagsList=[{}]
        simFlagsList[0]['storeQuantities']=[archive_pod_res]
    ns = bu.NumericalSolution.NS_base(bu.so,[bu.physics],[bu.numerics],bu.so.sList,bu.opts,simFlagsList=simFlagsList)
    return ns

def get_burgers_ns_red(name,T=0.1,nDTout=10,archive_pod_res=None,use_pod=True,use_hyper=False):
    import burgers_init as bu
    bu.physics.name=name
    bu.so.name = bu.physics.name
    #adjust default end time and number of output steps
    bu.T=T
    bu.nDTout=nDTout
    bu.DT=bu.T/float(bu.nDTout)
    bu.so.tnList = [i*bu.DT for i in range(bu.nDTout+1)]

    bu.numerics.SVD_basis_file='SVD_basis_truncated'
    bu.numerics.hyper_SVD_basis_file='Fn_SVD_basis_truncated'
    bu.numerics.hyper_indices_file = 'DEIM_indices'
    bu.numerics.hyper_Q_file = 'Q_DEIM_truncated'

    bu.numerics.use_hyper = use_hyper
    #use_hyper is dummy by default if use_pod == True
    if not use_pod:
	bu.numerics.multilevelNonlinearSolver = bu.NonlinearSolvers.POD_HyperReduced_Newton
	bu.numerics.levelNonlinearSolver = bu.NonlinearSolvers.POD_HyperReduced_Newton
    else:
	bu.numerics.multilevelNonlinearSolver = bu.NonlinearSolvers.POD_Newton
	bu.numerics.levelNonlinearSolver = bu.NonlinearSolvers.POD_Newton
    #request archiving of spatial residuals ...
    bu.numerics.tolFac = 0.0#relative tolerance
    bu.numerics.nl_atol_res = 1.0e-4 #nonlinear solver rtolerance
    simFlagsList=None
    if archive_pod_res is not None:
        simFlagsList=[{}]
        simFlagsList[0]['storeQuantities']=[archive_pod_res]
    ns = bu.NumericalSolution.NS_base(bu.so,[bu.physics],[bu.numerics],bu.so.sList,bu.opts,simFlagsList=simFlagsList)
    return ns

def test_burgers_run():
    """
    Aron's favority smoke test to see if burgers runs
    """
    ns = get_burgers_ns("test_run",T=0.1,nDTout=10)
    failed = ns.calculateSolution("run_smoke_test")
    assert not failed

def test_nonlin_residual_split():
    """
    just tests that R=R_l+R_n for random dof vector in [0,1]

    Here R_l and R_n are the linear and nonlinear portions of the residual
    """
    ns = get_burgers_ns("test_res_split",T=0.05,nDTout=5)
    failed = ns.calculateSolution("run_res_test")
    assert not failed
    #storage for residuals
    model = ns.modelList[0].levelModelList[-1]
    u_tmp = np.random.random(model.u[0].dof.shape)
    res = np.zeros(model.u[0].dof.shape,'d')
    res_n = res.copy(); res_l=res.copy()
    model.getResidual(u_tmp,res)
    model.getNonlinearResidual(u_tmp,res_n)
    model.getLinearResidual(u_tmp,res_l)
    
    res_n += res_l
    npt.assert_almost_equal(res,res_n)

def test_nonlin_jacobian_split():
    """
    just tests that J=J_l+J_n for random dof vector in [0,1]

    Here J_l and J_n are the linear and nonlinear portions of the Jacobian
    """
    ns = get_burgers_ns("test_jac_split",T=0.1,nDTout=10)
    failed = ns.calculateSolution("run_jac_test")
    assert not failed
    #storage for residuals
    model = ns.modelList[0].levelModelList[-1]
    u_tmp = np.random.random(model.u[0].dof.shape)
    res = np.zeros(model.u[0].dof.shape,'d')
    res_n = res.copy(); res_l=res.copy()
    model.getResidual(u_tmp,res)
    model.getNonlinearResidual(u_tmp,res_n)
    model.getLinearResidual(u_tmp,res_l)
    
    res_n += res_l
    npt.assert_almost_equal(res,res_n)
    #sparse storage for nonlinear and linear jacobians
    nonlin_J = model.initializeNonlinearJacobian()
    nonlin_J_rowptr,nonlin_J_colind,nonlin_J_nzval = nonlin_J.getCSRrepresentation()

    npt.assert_almost_equal(nonlin_J_rowptr,model.rowptr)
    npt.assert_almost_equal(nonlin_J_colind,model.colind)
    
    lin_J = model.initializeLinearJacobian()
    lin_J_rowptr,lin_J_colind,lin_J_nzval = lin_J.getCSRrepresentation()

    npt.assert_almost_equal(lin_J_rowptr,model.rowptr)
    npt.assert_almost_equal(lin_J_colind,model.colind)

    model.getLinearJacobian(lin_J)
    model.getNonlinearJacobian(nonlin_J)
    model.getJacobian(model.jacobian)

    nonlin_J_nzval += lin_J_nzval
    npt.assert_almost_equal(nonlin_J_nzval,model.nzval)
    
    
def test_res_archive():
    """
    smoke test if numerical solution can archive nonlinear residuals' to xdmf  
    """
    ns = get_burgers_ns("test_nonlin_res_archive",T=0.1,nDTout=10,archive_pod_res='pod_residuals_linnonlin')
    
    nonlin_failed = ns.calculateSolution("run_nonlin_res_smoke_test")
    assert not nonlin_failed


def test_svd_nonlin_res(file_prefix='Fn'):
    """
    test SVD decomposition of spatial residuals by generating SVD, saving to file and reloading
    """
    from proteus.deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_nonlinear_res",T=0.1,nDTout=10,archive_pod_res='pod_residuals_linnonlin')
    
    failed = ns.calculateSolution("run_svd_nonlinear_res")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_nonlinear_res",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'nonlinear_residual0',file_prefix)

    S_svd = np.dot(U,np.dot(np.diag(s),V))
    #now load back in and test
    S = read_snapshots(archive,len(ns.tnList),'nonlinear_residual0')

    npt.assert_almost_equal(S,S_svd)

def test_svd_linear_res(file_prefix='Fl'):
    """
    test SVD decomposition of mass residuals by generating SVD, saving to file and reloading
    """
    from proteus.deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_linear_res",T=0.1,nDTout=10,archive_pod_res='pod_residuals_linnonlin')
    
    failed = ns.calculateSolution("run_svd_linear_res")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_linear_res",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'linear_residual0',file_prefix)

    S_svd = np.dot(U,np.dot(np.diag(s),V))
    #now load back in and test
    S = read_snapshots(archive,len(ns.tnList),'linear_residual0')

    npt.assert_almost_equal(S,S_svd)

def test_svd_sol(file_prefix='sol'):
    """
    test SVD decomposition of solution by generating SVD, saving to file and reloading
    """
    from proteus.deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_sol",T=0.1,nDTout=10,archive_pod_res='pod_residuals_nonlin')
    
    failed = ns.calculateSolution("run_svd_sol")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_sol",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'u',file_prefix)

    S_svd = np.dot(U,np.dot(np.diag(s),V))
    #now load back in and test
    S = read_snapshots(archive,len(ns.tnList),'u')

    npt.assert_almost_equal(S,S_svd)

def test_deim_indices():
    """
    Taking a basis generated from snapshots
    and tests that get the deim algorithm returns all of the indices 
    """
    import os
    basis_file='Fn_SVD_basis'
    if not os.path.isfile(basis_file):
        test_svd_nonlin_res(file_prefix='Fn')
    U = np.loadtxt(basis_file)
    assert U.shape[1] > 1
    from proteus.deim_utils import calculate_deim_indices
    rho_half = calculate_deim_indices(U[:,:U.shape[1]/2])
    assert rho_half.shape[0] == U.shape[1]/2

    rho = calculate_deim_indices(U)    
    assert rho.shape[0] == U.shape[1]

    rho_uni = np.unique(rho)
    assert rho_uni.shape[0] == rho.shape[0]

    for i in range(rho_half.shape[0]):
	assert rho_half[i] == rho[i]

def deim_approx(T=0.1,nDTout=10,m=5,m_linear=5):
    """
    Follow basic setup for DEIM approximation
    - generate a burgers solution, saving nonlinear and linear residuals
    - generate SVDs for snapshots
    - for both residuals Fn and Fl
      - pick $m$, dimension for snapshot reduced basis $\mathbf{U}_m$
      - call DEIM algorithm to determine $\vec \rho$ and compute projection matrix 
        $\mathbf{P}_F=\mathbf{U}_m(\mathbf{P}^T\mathbf{U}_m)^{-1}$

    For selected timesteps
    - extract fine grid solution from archive, $\vec y$
    - for both residuals F=F_s and F_m
      - evaluate $\vec F(\vec y)$ at indices in $\vec \rho \rightarrow \vec c$
      - apply DEIM interpolant $\tilde{\vec F} = \mathbf{P}_F\vec c$
      - compute error $\|F-\tilde{\vec F}\|
    - visualize
    """

    from proteus.deim_utils import read_snapshots,generate_svd_decomposition
    ##run fine grid problem
    ns = get_burgers_ns("test_deim_approx",T=T,nDTout=nDTout,archive_pod_res='pod_residuals_linnonlin')
    
    failed = ns.calculateSolution("run_deim_approx")
    assert not failed

    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_deim_approx",readOnly=True)
    ##perform SVD on spatial residual
    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'nonlinear_residual0','Fn')
    U_l,s_l,V_l=generate_svd_decomposition(archive,len(ns.tnList),'linear_residual0','Fl')

    from proteus.deim_utils import deim_alg
    ##calculate DEIM indices and projection matrix
    rho,PF = deim_alg(U,m)
    #also 'linear' term
    rho_l,PF_l = deim_alg(U_l,m_linear)

    ##for comparison, grab snapshots of solution and residual
    Su = read_snapshots(archive,len(ns.tnList),'u')
    Sn = read_snapshots(archive,len(ns.tnList),'nonlinear_residual0')
    Sl = read_snapshots(archive,len(ns.tnList),'linear_residual0')

    steps_to_test = np.arange(len(ns.tnList)) 
    errors = np.zeros(len(steps_to_test),'d'); errors_lin = errors.copy()

    F_deim = np.zeros((Sn.shape[0],len(steps_to_test)),'d')
    Fl_deim = np.zeros((Sl.shape[0],len(steps_to_test)),'d')
    for i,istep in enumerate(steps_to_test):
        #solution on the fine grid
        u = Su[:,istep]
        #nonlinear residual evaluated from fine grid
        F = Sn[:,istep]
        #deim approximation on the fine grid 
        F_deim[:,istep] = np.dot(PF,F[rho])
        errors[i] = np.linalg.norm(F-F_deim[:,istep])
        #repeat for linear residual
        Fl= Sl[:,istep]
        #deim approximation on the fine grid 
        Fl_deim[:,istep] = np.dot(PF_l,Fl[rho_l])
        errors_lin[i] = np.linalg.norm(Fl-Fl_deim[:,istep])
    #
    np.savetxt("deim_approx_errors_nonlinear_test_T={0}_nDT={1}_m={2}.dat".format(T,nDTout,m),errors)
    np.savetxt("deim_approx_errors_linear_test_T={0}_nDT={1}_m={2}.dat".format(T,nDTout,m_linear),errors_lin)
        
    return errors,errors_lin,F_deim,Fl_deim

def test_deim_approx_full(tol=1.0e-12):
    """
    check that get very small error if use full basis 
    """
    T = 0.1; nDTout=10; m=nDTout+1
    errors,errors_linear,F_deim,Fl_deim = deim_approx(T=T,nDTout=nDTout,m=m,m_linear=m)
    assert errors.min() < tol
    assert errors_linear.min() < tol

def pod_hyper_pod_run(T=0.1,nDTout=10,m_sol=5,m=5):
    """
    Follow basic setup for DEIM approximation
    - generate a burgers solution, saving solution and nonlinear residuals
    - generate SVDs for snapshots
    - for residuals Fn
      - pick $m$, dimension for snapshot reduced basis $\mathbf{U}_m$
      - call DEIM algorithm to determine $\vec \rho$ and compute projection matrix 
        $\mathbf{P}_F=\mathbf{U}_m(\mathbf{P}^T\mathbf{U}_m)^{-1}$

    For all timesteps
    - extract both reduced solutions from archive
    - calculate the norms of the error between the two
    """

    from proteus.deim_utils import read_snapshots,generate_svd_decomposition
    ##run fine grid problem
    ns = get_burgers_ns("test_pod_hyper_pod_fine",T=T,nDTout=nDTout,archive_pod_res='pod_residuals_linnonlin')
    
    failed = ns.calculateSolution("run_pod_hyper_pod_fine")
    assert not failed

    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_pod_hyper_pod_fine",readOnly=True)
    ##perform SVD on spatial residual
    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'u','sol')
    Un,sn,Vn=generate_svd_decomposition(archive,len(ns.tnList),'nonlinear_residual0','Fn')

    from proteus.deim_utils import deim_alg
    ##calculate DEIM indices and projection matrix
    rho,PF = deim_alg(Un,m)

    np.savetxt('SVD_basis_truncated', U[:,0:m_sol], delimiter=' ')
    np.savetxt('Fn_SVD_basis_truncated', Un[:,0:m], delimiter=' ')
    np.savetxt('DEIM_indices', rho, delimiter=' ',fmt='%d')
    np.savetxt('Q_DEIM_truncated', PF, delimiter=' ')

    ##reduced order models below
    ##pure POD first
    ns_red = get_burgers_ns_red("test_pod_hyper_pod_red",T=T,nDTout=nDTout,archive_pod_res=None,use_pod=True,use_hyper=False)

    failed = ns_red.calculateSolution("run_pod_hyper_pod_red")
    assert not failed

    archive_red = Archiver.XdmfArchive(".","test_pod_hyper_pod_red",readOnly=True)
    ##for comparison, grab snapshots of solution
    Su = read_snapshots(archive_red,len(ns_red.tnList),'u')

    ##POD-DEIM second
    ns_red1 = get_burgers_ns_red("test_pod_hyper_pod_red1",T=T,nDTout=nDTout,archive_pod_res=None,use_pod=False,use_hyper=False)

    failed = ns_red1.calculateSolution("run_pod_hyper_pod_red1")
    assert not failed

    archive_red1 = Archiver.XdmfArchive(".","test_pod_hyper_pod_red1",readOnly=True)
    ##for comparison, grab snapshots of solution
    Su1 = read_snapshots(archive_red1,len(ns_red1.tnList),'u')

    errors = np.zeros(Su.shape[0],'d');
    for i in range(nDTout):
	errors[i] = np.linalg.norm(Su[:,i]-Su1[:,i])
    #
    np.savetxt("pod_poddeim_errors_T={0}_nDT={1}_msol={2}_m={3}.dat".format(T,nDTout,m_sol,m),errors)
        
    return errors

def test_pod_hyper_pod(tol=1.0e-14):
    """
    check that pod-deim = pod is use_hyper = False 
    """
    T = 0.1; nDTout=10; m_sol=5; m=5
    errors = pod_hyper_pod_run(T=T,nDTout=nDTout,m_sol=m_sol,m=m)
    assert errors.min() < tol

def deim_run():
    """
    Follow basic setup for DEIM approximation
    - generate a burgers solution, saving solution and nonlinear residuals
    - generate SVDs for snapshots
    - for residuals Fn
      - pick $m=10$, dimension for snapshot reduced basis $\mathbf{U}_m$
      - call DEIM algorithm to determine $\vec \rho$ and compute projection matrix 
        $\mathbf{P}_F=\mathbf{U}_m(\mathbf{P}^T\mathbf{U}_m)^{-1}$
    """

    from proteus.deim_utils import read_snapshots,generate_svd_decomposition
    ##run fine grid problem
    ##these four values are rigidly fixed!!!!!
    T = 1; nDTout = 100; m_sol = 10; m = 10
    ns = get_burgers_ns("test_fine_run",T=T,nDTout=nDTout,archive_pod_res='pod_residuals_linnonlin')
    
    failed = ns.calculateSolution("run_fine_run")
    assert not failed

    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_fine_run",readOnly=True)
    ##perform SVD on nonlinear residual
    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'u','sol')
    Un,sn,Vn=generate_svd_decomposition(archive,len(ns.tnList),'nonlinear_residual0','Fn')

    from proteus.deim_utils import deim_alg
    ##calculate DEIM indices and projection matrix
    rho,PF = deim_alg(Un,m)

    np.savetxt('SVD_basis_truncated', U[:,0:m_sol], delimiter=' ')
    np.savetxt('Fn_SVD_basis_truncated', Un[:,0:m], delimiter=' ')
    np.savetxt('DEIM_indices', rho, delimiter=' ',fmt='%d')
    np.savetxt('Q_DEIM_truncated', PF, delimiter=' ')

    ##reduced order models below
    ##pure POD first
    ns_red = get_burgers_ns_red("test_deim_run",T=T,nDTout=nDTout,archive_pod_res=None,use_pod=False,use_hyper=True)

    failed = ns_red.calculateSolution("run_deim_run")
    assert not failed

    archive_red = Archiver.XdmfArchive(".","test_deim_run",readOnly=True)
    ##for comparison, grab snapshots of solution
    Su = read_snapshots(archive_red,len(ns_red.tnList),'u')
    np.savetxt("proteus_solution_T={0}_nDT={1}_msol={2}_m={3}.dat".format(T,nDTout,m_sol,m),Su)

    return

def test_deim_impl():
    """
    compare output from matlab POD-DEIM with proteus POD-DEIM
    """
    T = 1; nDTout=100; m_sol=10; m=10
    #matlab solution: m_sol = 10, m = 10
    solm = np.loadtxt("S_matlab")
    solm = solm[:,-1]
    #proteus solution
    deim_run()
    sol = np.loadtxt("proteus_solution_T=1_nDT=100_msol=10_m=10.dat")#.format(T,nDTout,m_sol,m))
    sol = sol[:,-1]

    diff = solm-sol
    imax = np.argmax(np.absolute(diff))
    assert np.absolute(diff.flat[imax]) < 1.0e-10, "got larger value than expected at index {0}: {1}".format(imax,diff.flat[imax])
    npt.assert_almost_equal(sol,solm)
    
if __name__ == "__main__":
    from proteus import Comm
    comm = Comm.init()
    import nose
    #nose.main(defaultTest='test_deim:test_burgers_run')
    nose.main(defaultTest='test_deim:test_nonlin_residual_split')
