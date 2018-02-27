from nose.tools import eq_ as eq
from nose.tools import ok_ as ok
import subprocess
import os
import pytest

@pytest.mark.slowTest
def test_workflowPUMI(verbose=0):
    """Test serial workflow: load model and mesh, solve within proteus 
    estimate error, adapt, solve again. It's not so important if the
    problem isn't setup properly"""
    subprocess.call("parun couette_so.py",shell=True )
    assert(True)
if __name__ == '__main__':
    import nose
    nose.main(defaultTest='test_workflowPUMI:test_workflowPUMI')


