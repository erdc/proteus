#!/usr/bin/env python

"""Detect MATLAB on the system, and configure both Proteus and HashStack to use it

This is designed to work on Linux, OS X, and Cygwin.  Native Windows is currently not supported.

This relies on the marked_yaml from a HashDist installation being available
"""

import sys
import os.path
import subprocess
try:
    import yaml
except ImportError:
    # assume hashdist source is available in this repository
    sys.path.append('./hashdist')
    import hashdist.deps.yaml as yaml

def detect_matlab(where=''):
    """Return path to MATLAB binary.

    First look in where, then on the PATH using which"""

    if where:
        for loc in (where, os.path.join(where, 'matlab')):
            abs_loc = os.path.abspath(loc)
            if os.path.isfile(abs_loc):
                return abs_loc
        else:
            raise Exception("Unable to locate MATLAB in {}".format(where))

    # No where specified, try the PATH

    loc = subprocess.check_output(['which', 'matlab'])
    return loc


def set_stack_parameters(profile_file, path_to_matlab, package_dict):
    """ Given a HashStack profile, set the MATLAB path and append the listed packages
    """

    with open(profile_file) as profile_yaml:
        profile = yaml.safe_load(profile_yaml)

    sys.stdout.write('Appending {} to packages\n'.format(package_dict))
    profile['packages'].update(package_dict)
    sys.stdout.write('Appending HOST_MATLAB = {} to parameters\n'.format(path_to_matlab))
    if 'parameters' in profile:
        if profile['parameters'] is None:
            profile['parameters'] = {}
        profile['parameters']['HOST_MATLAB'] = path_to_matlab
    else:
        profile['parameters'] = {'HOST_MATLAB':path_to_matlab}

    with open(profile_file, 'w') as profile_yaml:
        yaml.safe_dump(profile, profile_yaml)

def setup_proteus(profile_file, path_to_matlab=''):
    """ Detects MATLAB on PATH, sets HOST_MATLAB, and adds pymatbridge to packages in profile_file
    """

    path_to_matlab = detect_matlab(path_to_matlab)
    package_dict = {'pymatbridge':None,'matlab':{'use':'host-matlab'}}
    set_stack_parameters(profile_file, path_to_matlab, package_dict)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.stderr.write("setupmatlab.py: You must specify a profile file to modify\n")
        sys.exit(1)
    profile_file = sys.argv[1]
    try:
        if len(sys.argv) > 2:
            setup_proteus(profile_file, path_to_matlab=sys.argv[2])
        else:
            setup_proteus(profile_file)
    except subprocess.CalledProcessError:
        sys.stderr.write('setupmatlab.py: Unable to find MATLAB using "which matlab", please add it to your PATH\n')
        sys.exit(1)
    except Exception:
        sys.stderr.write('setupmatlab.py: Unable to find MATLAB in the given path: {}\n'.format(sys.argv[2]))
        sys.exit(1)
