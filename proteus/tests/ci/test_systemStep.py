#!/usr/bin/env python
"""
Test module for fixed time stepping (serial)

This module solves equations of the form

.. _math::

  u_t + \nabla \cdot \left( u \mathbf{v} - a(x) \nabla u \right) = 0

"""
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
Profiling.logLevel=7
Profiling.verbose=False
from petsc4py import PETSc
import numpy as np
import numpy.testing as npt
import ladr_2d_p
import ladr_2d_n


def test_minModelStep_stepExactTrue():
    reload(ladr_2d_p)
    reload(ladr_2d_n)
    pList = [ladr_2d_p]
    nList = [ladr_2d_n]
    reload(default_so)
    so = default_so
    so.name = pList[0].name = "ladr"
    so.tnList = nList[0].tnList
    so.systemStepControllerType = SplitOperator.Sequential_MinModelStep
    so.systemStepExact=True
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].runCFL=0.33
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('ladr_minModelStep_stepExactTrue')
    assert ns.tCount + 1 == len(so.tnList), "wrong number of archvie steps "+`ns.tCount`
    assert ns.modelList[0].solver.solverList[0].solveCalls == 40, "wrong number of steps "+`ns.modelList[0].solver.solverList[0].solveCalls`
    archiveTimes=[]
    for t in ns.ar[0].treeGlobal.iter('Time'):
        archiveTimes.append(t.attrib['Value'])
    archiveTimesCorrect = so.tnList
    npt.assert_almost_equal(np.array(archiveTimes,'d'), np.array(archiveTimesCorrect,'d'))
    del ns

def test_minModelStep_stepExactFalse():
    reload(ladr_2d_p)
    reload(ladr_2d_n)
    pList = [ladr_2d_p]
    nList = [ladr_2d_n]
    reload(default_so)
    so = default_so
    so.name = pList[0].name = "ladr"
    so.tnList = nList[0].tnList
    so.systemStepControllerType = SplitOperator.Sequential_MinModelStep
    so.systemStepExact=False
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].runCFL=0.33
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('ladr_minModelStep_stepExactFalse')
    assert ns.tCount + 1 == len(so.tnList), "wrong number of archvie steps " +`ns.tCount`
    assert ns.modelList[0].solver.solverList[0].solveCalls == 34, "wrong number of steps "+`ns.modelList[0].solver.solverList[0].solveCalls`
    archiveTimes=[]
    for t in ns.ar[0].treeGlobal.iter('Time'):
        archiveTimes.append(t.attrib['Value'])
    archiveTimesCorrect = ['0.0', '0.029516097303', '0.0516531702802', '0.0811692675832', '0.10330634056', '0.125443413538', '0.154959510841', '0.177096583818', '0.206612681121', '0.228749754098', '0.250886827075']
    npt.assert_almost_equal(np.array(archiveTimes,'d'), np.array(archiveTimesCorrect,'d'))
    del ns

def test_fixedStep_stepExactFalse():
    reload(ladr_2d_p)
    reload(ladr_2d_n)
    pList = [ladr_2d_p]
    nList = [ladr_2d_n]
    reload(default_so)
    so = default_so
    so.name = pList[0].name = "ladr"
    so.tnList = nList[0].tnList
    so.systemStepControllerType = SplitOperator.Sequential_FixedStep
    so.systemStepExact=False
    so.dt_system_fixed = 0.01
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].runCFL=0.33
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('ladr_minModelStep_stepExactFalse')
    assert ns.tCount + 1 == len(so.tnList), "wrong number of archvie steps " +`ns.tCount`
    assert ns.modelList[0].solver.solverList[0].solveCalls == 25, "wrong number of steps "+`ns.modelList[0].solver.solverList[0].solveCalls`
    archiveTimes=[]
    for t in ns.ar[0].treeGlobal.iter('Time'):
        archiveTimes.append(t.attrib['Value'])
    archiveTimesCorrect = []
    interval=0
    step=0
    t=so.tnList[0]
    archiveTimesCorrect.append(t)
    while interval < len(so.tnList)-1:
        t+= so.dt_system_fixed
        if t >= so.tnList[interval+1]:
            archiveTimesCorrect.append(t)
            interval+=1
    npt.assert_almost_equal(np.array(archiveTimes,'d'), np.array(archiveTimesCorrect,'d'))
    del ns

def test_fixedStep_stepExactTrue():
    reload(ladr_2d_p)
    reload(ladr_2d_n)
    pList = [ladr_2d_p]
    nList = [ladr_2d_n]
    reload(default_so)
    so = default_so
    so.name = pList[0].name = "ladr"
    so.tnList = nList[0].tnList
    so.systemStepControllerType = SplitOperator.Sequential_FixedStep
    so.systemStepExact=True
    so.dt_system_fixed = 0.01
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].runCFL=0.33
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('ladr_minModelStep_stepExactFalse')
    assert ns.tCount + 1 == len(so.tnList), "wrong number of archvie steps " +`ns.tCount`
    assert ns.modelList[0].solver.solverList[0].solveCalls == 30, "wrong number of steps "+`ns.modelList[0].solver.solverList[0].solveCalls`
    archiveTimes=[]
    for t in ns.ar[0].treeGlobal.iter('Time'):
        archiveTimes.append(t.attrib['Value'])
    archiveTimesCorrect = so.tnList
    npt.assert_almost_equal(np.array(archiveTimes,'d'), np.array(archiveTimesCorrect,'d'))
    del ns

def test_fixedStep_stepSimple():
    reload(ladr_2d_p)
    reload(ladr_2d_n)
    pList = [ladr_2d_p]
    nList = [ladr_2d_n]
    reload(default_so)
    so = default_so
    so.name = pList[0].name = "ladr"
    so.tnList = nList[0].tnList
    so.systemStepControllerType = SplitOperator.Sequential_FixedStep_Simple
    so.systemStepExact=True#should be ignored
    so.dt_system_fixed = 0.01#should be ignored
    so.sList=[default_s]
    opts.logLevel=7
    opts.verbose=True
    opts.profile=True
    opts.gatherArchive=True
    nList[0].runCFL=0.33
    nList[0].linearSolver=default_n.LU
    nList[0].multilevelLinearSolver=default_n.LU
    ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
    ns.calculateSolution('ladr_minModelStep_stepExactFalse')
    assert ns.tCount + 1 == len(so.tnList), "wrong number of archvie steps " +`ns.tCount`
    assert ns.modelList[0].solver.solverList[0].solveCalls == len(so.tnList)-1, "wrong number of steps "+`ns.modelList[0].solver.solverList[0].solveCalls`
    archiveTimes=[]
    for t in ns.ar[0].treeGlobal.iter('Time'):
        archiveTimes.append(t.attrib['Value'])
    archiveTimesCorrect = so.tnList
    npt.assert_almost_equal(np.array(archiveTimes,'d'), np.array(archiveTimesCorrect,'d'))
    del ns

if __name__ == '__main__':
    test_minModelStep_stepExactTrue() 
    test_minModelStep_stepExactFalse()
    test_fixedStep_stepExactFalse()
    test_fixedStep_stepExactTrue()
    Profiling.logEvent("Closing Log")
    try:
        Profiling.closeLog()
    except:
        pass
