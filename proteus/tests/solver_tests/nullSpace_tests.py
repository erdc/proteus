#!usr/bin/env python
""" Test module for null space """

from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel = 7
Profiling.verbose  = True
from nose.tools import ok_ as ok
from nose.tools import eq_ as eq
from petsc4py import PETSc as p4pyPETSc

# TODO - write test classes with setup and teardown

class stokesCavity_test():
    def setUp(self):
        pass
    def tearDown(self):
        pass
    

def stokesCavity_test():
    """ This function tests the driven cavity """
    import stokesDrivenCavity_2d_p
    import stokesDrivenCavity_2d_n
    pList = [stokesDrivenCavity_2d_p]
    nList = [stokesDrivenCavity_2d_n]
    so = default_so
    so.tnList = [0.,1.]
    so.name = pList[0].name
    so.sList = [default_s]
    opts.verbose = True
    opts.profile = True
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('test_NullSpace')
    


if __name__=='__main__':
    from proteus import Comm
    comm = Comm.init()
    stokesCavity_test()
