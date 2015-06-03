
def test_povgen():
    import difflib
    import urllib
    import subprocess
    import os
    urllib.urlretrieve(
        'http://dl.dropboxusercontent.com/u/26353144/floating_bar0.h5',
        'floating_bar0.h5')
    subprocess.check_call(['../../../scripts/povgen.py',
                           'floating_bar',
                           '-s',
                           '3'])
    # prefix='phi_0.000000_
    povfiles = []
    for i in range(3):
        filename = 'phi_0.000000_{0:04d}.pov'.format(i)
        with open(filename, 'r') as f:
            povfiles.append(f.readlines())
    urllib.urlretrieve(
        'http://dl.dropboxusercontent.com/u/26353144/phi_0.000000_000.tgz',
        'phi_0.000000_000.tgz')
    subprocess.check_call(['tar', 'xzf', 'phi_0.000000_000.tgz'])
    saved_povfiles = []
    for i in range(3):
        filename = 'phi_0.000000_{0:04d}.pov'.format(i)
        with open(filename, 'r') as f:
            saved_povfiles.append(f.readlines())
        assert saved_povfiles[i] == povfiles[i], \
            ''.join(list(difflib.unified_diff(saved_povfiles[i],
                                              povfiles[i],
                                              "archived",
                                              "test")))
    os.remove('phi_0.000000_000.tgz')
    os.remove('floating_bar0.h5')
