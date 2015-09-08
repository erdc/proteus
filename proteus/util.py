""" basic utilities in Proteus
"""

import os

def get_include_dir():
    return os.path.dirname(os.path.realpath(__file__))
