import os
import pytest

class TestIsosurface(object):

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
    
    @pytest.mark.skip
    def test_povgen(self):
        import difflib
        import subprocess
        import os
        subprocess.check_call(['povgen.py',
                               os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            'comparison_files',
                                            'floating_bar'),
                               '-s',
                               '3'])
        povfiles = []
        for i in range(3):
            filename = 'phi_t_0.000000_{0:04d}.pov'.format(i)
            with open(filename, 'r') as f:
                povfiles.append(f.readlines())
        subprocess.check_call(['tar', 'xzf', os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                          'comparison_files',
                                                          'phi_t_0.000000_000.tgz')])
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
