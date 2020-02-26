import subprocess
import os
import pytest

@pytest.mark.slowTest
def test_workflowPUMI(verbose=0):
    """Test serial workflow: load model and mesh, solve within proteus 
    estimate error, adapt, solve again. It's not so important if the
    problem isn't setup properly"""
    currentPath = os.path.dirname(os.path.abspath(__file__))
    runCommand = "cd "+currentPath+"; parun -l5 couette_so.py;"
    subprocess.check_call(runCommand,shell=True )
if __name__ == '__main__':
    import nose
    nose.main(defaultTest='test_workflowPUMI:test_workflowPUMI')


