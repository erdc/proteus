#!/usr/bin/env python

import subprocess
import os
import pytest

@pytest.mark.slowTest

def test_MeshAdaptRestart_generateMesh(verbose=0):
    """Generate Mesh for PUMI"""
    currentPath = os.path.dirname(os.path.abspath(__file__))
    runCommand = "cd "+currentPath+"; parun -C \"gen_mesh=True usePUMI=True adapt=0 fixedTimeStep=False\" -D \"baseline\" dambreak_Colagrossi_so.py;"
    subprocess.check_call(runCommand,shell=True)

def test_MeshAdaptRestart_adaptiveTime_BackwardEuler_baseline_withRedistancing(verbose=0):
    """Get baseline data"""
    currentPath = os.path.dirname(os.path.abspath(__file__))
    runCommand = "cd "+currentPath+"; parun -C \"gen_mesh=False usePUMI=True adapt=0 fixedTimeStep=False\" -D \"baseline\" dambreak_Colagrossi_so.py;"
    subprocess.check_call(runCommand,shell=True)

def test_MeshAdaptRestart_adaptiveTime_BackwardEuler_withRedistancing(verbose=0):
    """Test restart workflow"""
    currentPath = os.path.dirname(os.path.abspath(__file__))
    runCommand = "cd "+currentPath+"; parun -C \"gen_mesh=False usePUMI=True adapt=1 fixedTimeStep=False\" -D \"adapt_0\" dambreak_Colagrossi_so.py;"
    subprocess.check_call(runCommand,shell=True)

def test_MeshAdaptRestart_withRedistancing(verbose=0):
    #subprocess.call("diff normal adapt_0/pressure.csv baseline/pressure.csv")
    currentPath = os.path.dirname(os.path.abspath(__file__))
    with open(currentPath+'/baseline/pressure.csv') as file1, open(currentPath+'/adapt_0/pressure.csv') as file2:
        for line1, line2 in zip(file1, file2):
            if(line1 != line2):
                pytest.fail("pressure gauge values are not the same!\n")
            