#!/usr/bin/env python

"""
some tests for initial DEIM implementation in proteus
"""
import numpy as np
import numpy.testing as npt
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq

def get_burgers_ns(name,T=0.1,nDTout=10):
    import burgers_init as bu
    bu.physics.name=name
    bu.so.name = bu.physics.name
    #adjust default end time and number of output steps
    bu.T=T
    bu.nDTout=nDTout
    bu.DT=bu.T/float(bu.nDTout)
    bu.so.tnList = [i*bu.DT for i in range(bu.nDTout+1)]
    ns = bu.NumericalSolution.NS_base(bu.so,[bu.physics],[bu.numerics],bu.so.sList,bu.opts)
    return ns

def test_burger_run():
    ns = get_burgers_ns("test_run",T=0.1,nDTout=10)
    failed = ns.calculateSolution("run_smoke_test")
    assert not failed

def test_residual_split():
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

if __name__ == "__main__":
    from proteus import Comm
    comm = Comm.init()
    import nose
    nose.main(defaultTest='test_DEIM:test_residual_split')
