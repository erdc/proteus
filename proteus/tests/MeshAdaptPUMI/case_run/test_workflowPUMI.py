import subprocess
import os
import pytest

@pytest.mark.slowTest
def test_workflowPUMI(verbose=0):
    """Test serial workflow: load model and mesh, solve within proteus 
    estimate error, adapt, solve again. It's not so important if the
    problem isn't setup properly"""
    subprocess.check_call("parun -l5 couette_so.py",shell=True )
if __name__ == '__main__':
    import nose
    nose.main(defaultTest='test_workflowPUMI:test_workflowPUMI')


