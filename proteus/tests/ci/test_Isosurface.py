import os

class TestIsosurface():

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        pass

    def teardown_method(self,method):
        [os.remove(f) for f in os.listdir(".") if f.endswith(".pov")]
        FileList = ['isostats',
                    'proteus_default.log',
                    'proteus.inc',
                    'proteus.log']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
        
    def test_povgen(self):
        import difflib
        import urllib
        import subprocess
        import os
        urllib.urlretrieve(
            'https://dl.dropbox.com/s/tjkj3ella3ntv75/floating_bar.h5',
            'floating_bar.h5')
        subprocess.check_call(['povgen.py',
                               'floating_bar',
                               '-s',
                               '3'])
        povfiles = []
        for i in range(3):
            filename = 'phi_t_0.000000_{0:04d}.pov'.format(i)
            with open(filename, 'r') as f:
                povfiles.append(f.readlines())
        urllib.urlretrieve(
            'https://dl.dropbox.com/s/um635cllwdjxajc/phi_t_0.000000_000.tgz',
            'phi_t_0.000000_000.tgz')
        subprocess.check_call(['tar', 'xzf', 'phi_t_0.000000_000.tgz'])
        saved_povfiles = []
        for i in range(3):
            filename = 'phi_t_0.000000_{0:04d}.pov'.format(i)
            with open(filename, 'r') as f:
                saved_povfiles.append(f.readlines())
            assert saved_povfiles[i] == povfiles[i], \
                ''.join(list(difflib.unified_diff(saved_povfiles[i],
                                                  povfiles[i],
                                                  "archived",
                                                  "test")))
        os.remove('phi_t_0.000000_000.tgz')
        os.remove('floating_bar.h5')
