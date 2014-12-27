#!/usr/bin/env python

"""
some tests for initial DEIM implementation in proteus
"""
import numpy as np
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

def get_burgers_ns(name,T=0.1,nDTout=10,archive_space_res=False):
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
    if archive_space_res:
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
    ns = get_burgers_ns("test_space_res_archive",T=0.1,nDTout=10,archive_space_res=True)
    
    failed = ns.calculateSolution("run_space_res_smoke_test")
    assert not failed

def test_svd_space_res(file_prefix='F_s'):
    """
    test SVD decomposition of spatial residuals by generating SVD, saving to file and reloading
    """
    from deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_space_res",T=0.1,nDTout=10,archive_space_res=True)
    
    failed = ns.calculateSolution("run_svd_space_res")
    assert not failed
    from proteus import Archiver
    archive = Archiver.XdmfArchive(".","test_svd_space_res",readOnly=True)

    U,s,V=generate_svd_decomposition(archive,len(ns.tnList),'spatial_residual0','F_s')

    S_svd = np.dot(U,np.dot(np.diag(s),V))
    #now load back in and test
    S = read_snapshots(archive,len(ns.tnList),'spatial_residual0')

    npt.assert_almost_equal(S,S_svd)

def test_svd_soln():
    """
    test SVD decomposition of solution by generating SVD, saving to file and reloading
    """
    from deim_utils import read_snapshots,generate_svd_decomposition

    ns = get_burgers_ns("test_svd_soln",T=0.1,nDTout=10,archive_space_res=True)
    
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
    from deim_utils import calculate_deim_indices
    rho_half = calculate_deim_indices(U[:,:U.shape[1]/2])
    assert rho_half.shape[0] == U.shape[1]/2

    rho = calculate_deim_indices(U)    
    assert rho.shape[0] == U.shape[1]

    rho_uni = np.unique(rho)
    assert rho_uni.shape[0] == rho.shape[0]

if __name__ == "__main__":
    from proteus import Comm
    comm = Comm.init()
    import nose
    nose.main(defaultTest='test_deim:test_deim_indices')
