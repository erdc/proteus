#!/usr/bin/env python
"""
some tests for initial DEIM implementation in proteus
"""
import numpy as np
import numpy.testing as npt
import sys, os
from proteus import deim_utils

def get_burgers_ns(name,T=0.1,nDTout=10,archive_pod_res=False):
    sys.path.append(os.path.dirname(__file__))
    import burgers_init as bu
    sys.path.remove(os.path.dirname(__file__))
    bu.physics.name=name
    bu.so.name = bu.physics.name
    #adjust default end time and number of output steps
    bu.T=T
    bu.nDTout=nDTout
    bu.DT=bu.T/float(bu.nDTout)
    bu.so.tnList = [i*bu.DT for i in range(bu.nDTout+1)]
    #request archiving of spatial residuals ...
    simFlagsList=None
    if archive_pod_res:
        simFlagsList=[{}]
        simFlagsList[0]['storeQuantities']=['pod_residuals']
    ns = bu.NumericalSolution.NS_base(bu.so,[bu.physics],[bu.numerics],bu.so.sList,bu.opts,simFlagsList=simFlagsList)
    return ns

def test_burger_run():
    """
    Aron's favority smoke test to see if burgers runs
    """
    ns = get_burgers_ns("test_run",T=0.1,nDTout=10)
    failed = ns.calculateSolution("run_smoke_test")
    assert not failed

def test_residual_split():
    """
    just tests that R=R_s+R_t for random dof vector in [0,1]

    Here R_s and R_t are the spatial and mass residuals
    """
    ns = get_burgers_ns("test_res_split",T=0.05,nDTout=5)
    failed = ns.calculateSolution("run_res_test")
    assert not failed
    #storage for residuals
    model = ns.modelList[0].levelModelList[-1]
    u_tmp = np.random.random(model.u[0].dof.shape)
    res = np.zeros(model.u[0].dof.shape,'d')
    res_s = res.copy(); res_t=res.copy()
    model.getResidual(u_tmp,res)
    model.getSpatialResidual(u_tmp,res_s)
    model.getMassResidual(u_tmp,res_t)
    
    res_t += res_s
    npt.assert_almost_equal(res,res_t)

def test_res_archive():
    """
    smoke test if numerical solution can archive 'spatial residuals' to xdmf  
    """
    ns = get_burgers_ns("test_space_res_archive",T=0.1,nDTout=10,archive_pod_res=True)
    
    failed = ns.calculateSolution("run_space_res_smoke_test")
    assert not failed

def test_svd_space_res(file_prefix='F_s'):
    """
    test SVD decomposition of spatial residuals by generating SVD, saving to file and reloading
    """
    from proteus.deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_space_res",T=0.1,nDTout=10,archive_pod_res=True)
    
    failed = ns.calculateSolution("run_svd_space_res")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_space_res",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'spatial_residual0',file_prefix)

    S_svd = np.dot(U,np.dot(np.diag(s),V))
    #now load back in and test
    S = read_snapshots(archive,len(ns.tnList),'spatial_residual0')

    npt.assert_almost_equal(S,S_svd)

def test_svd_mass_res(file_prefix='F_m'):
    """
    test SVD decomposition of mass residuals by generating SVD, saving to file and reloading
    """
    from proteus.deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_mass_res",T=0.1,nDTout=10,archive_pod_res=True)
    
    failed = ns.calculateSolution("run_svd_mass_res")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_mass_res",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'mass_residual0',file_prefix)

    S_svd = np.dot(U,np.dot(np.diag(s),V))
    #now load back in and test
    S = read_snapshots(archive,len(ns.tnList),'mass_residual0')

    npt.assert_almost_equal(S,S_svd)

def test_svd_soln():
    """
    test SVD decomposition of solution by generating SVD, saving to file and reloading
    """
    from proteus.deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_soln",T=0.1,nDTout=10,archive_pod_res=True)
    
    failed = ns.calculateSolution("run_svd_soln")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_soln",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'u','soln')

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
    basis_file='F_s_SVD_basis'
    if not os.path.isfile(basis_file):
        test_svd_space_res(file_prefix='F_s')
    U = np.loadtxt(basis_file)
    from proteus.deim_utils import calculate_deim_indices
    rho_half = calculate_deim_indices(U[:,:U.shape[1]//2])
    assert rho_half.shape[0] == U.shape[1]//2

    rho = calculate_deim_indices(U)    
    assert rho.shape[0] == U.shape[1]

    rho_uni = np.unique(rho)
    assert rho_uni.shape[0] == rho.shape[0]

def deim_approx(T=0.1,nDTout=10,m=5,m_mass=5):
    """
    Follow basic setup for DEIM approximation
    - generate a burgers solution, saving spatial and 'mass' residuals
    - generate SVDs for snapshots
    - for both residuals F_s and F_m
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
    ns = get_burgers_ns("test_deim_approx",T=T,nDTout=nDTout,archive_pod_res=True)
    
    failed = ns.calculateSolution("run_deim_approx")
    assert not failed

    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_deim_approx",readOnly=True)
    ##perform SVD on spatial residual
    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'spatial_residual0','F_s')
    U_m,s_m,V_m=generate_svd_decomposition(archive,len(ns.tnList),'mass_residual0','F_m')

    from proteus.deim_utils import deim_alg
    ##calculate DEIM indices and projection matrix
    rho,PF = deim_alg(U,m)
    #also 'mass' term
    rho_m,PF_m = deim_alg(U_m,m_mass)

    ##for comparison, grab snapshots of solution and residual
    Su = read_snapshots(archive,len(ns.tnList),'u')
    Sf = read_snapshots(archive,len(ns.tnList),'spatial_residual0')
    Sm = read_snapshots(archive,len(ns.tnList),'mass_residual0')

    steps_to_test = np.arange(len(ns.tnList)) 
    errors = np.zeros(len(steps_to_test),'d'); errors_mass = errors.copy()

    F_deim = np.zeros((Sf.shape[0],len(steps_to_test)),'d')
    Fm_deim = np.zeros((Sf.shape[0],len(steps_to_test)),'d')
    for i,istep in enumerate(steps_to_test):
        #solution on the fine grid
        u = Su[:,istep]
        #spatial residual evaluated from fine grid
        F = Sf[:,istep]
        #deim approximation on the fine grid 
        F_deim[:,istep] = np.dot(PF,F[rho])
        errors[i] = np.linalg.norm(F-F_deim[:,istep])
        #repeat for 'mass residual'
        Fm= Sm[:,istep]
        #deim approximation on the fine grid 
        Fm_deim[:,istep] = np.dot(PF_m,Fm[rho])
        errors_mass[i] = np.linalg.norm(Fm-Fm_deim[:,istep])
    #
    np.savetxt("deim_approx_errors_space_test_T={0}_nDT={1}_m={2}.dat".format(T,nDTout,m),errors)
    np.savetxt("deim_approx_errors_mass_test_T={0}_nDT={1}_m={2}.dat".format(T,nDTout,m_mass),errors_mass)
        
    return errors,errors_mass,F_deim,Fm_deim

def test_deim_approx_full(tol=1.0e-12):
    """
    check that get very small error if use full basis 
    """
    T = 0.1; nDTout=10; m=nDTout+1
    errors,errors_mass,F_deim,Fm_deim = deim_approx(T=T,nDTout=nDTout,m=m,m_mass=m)
    assert errors.min() < tol
    assert errors_mass.min() < tol

