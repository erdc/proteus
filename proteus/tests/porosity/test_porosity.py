import subprocess
import os
import pytest

@pytest.mark.porosity

def teardown_method():
    """ Tear down function """
    FileList = ['sloshing.xmf',
                'sloshing.h5',
                'forceHistory_p.txt',
                'forceHistory_v.txt',
                'momentHistory.txt',
                'wettedAreaHistory.txt',
                'timeHistory.txt',
                'sloshing.log'
    ]
    for file in FileList:
        if os.path.isfile(file):
            os.remove(file)
        else:
            pass

def test_runTwoPhaseFlow_wPorosity_2D():
    teardown_method()
    """Test restart workflow"""
    currentPath = os.path.dirname(os.path.abspath(__file__))
    runCommand = "cd "+currentPath+"; parun --TwoPhaseFlow -l5 -v -C \"T=20. he=0.025 openTop=False genMesh=False\" sloshing.py;"
    subprocess.check_call(runCommand,shell=True)
    
def test_massConservation_wPorosity_2D(mcorr_nl_atol_res=1e-10):
    currentPath = os.path.dirname(os.path.abspath(__file__))
    file=open(currentPath+'/sloshing.log','r')
    lines=file.readlines()
    mcorr=[]
    nSteps=0
    for i,j in enumerate(lines):
        if "Phase 0 mass (consistent) after mass correction" in j:
            mcorr.append([float(lines[i+1].split()[7][:-1]),float(j.split()[-1])])
            nSteps +=1
    if (mcorr[0][1]+nSteps*mcorr_nl_atol_res)<=mcorr[-1][-1] or (mcorr[0][1]-nSteps*mcorr_nl_atol_res)>=mcorr[-1][-1]:
        pytest.fail("Mass is not being conserved within nonlinear residual tolerance for MCorr. Absolute mass error: {0:e}\n".format(abs(mcorr[-1][-1]-mcorr[0][1])))