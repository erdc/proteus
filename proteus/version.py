"""
Determine versions of Proteus/HashDist/HashStack by inspecting version files
"""

import os
from os.path import join as pjoin
import sys
import subprocess

git_cmd = ['git', 'log', '-1', '--pretty=%H']

# hashdist/hashstack are just files

hashdist = 'unknown'
hashstack = 'unknown'
proteus = 'unknown'

try:
    with open(pjoin(sys.prefix, 'hashdist_version.txt')) as f:
        hashdist = f.readline().strip()
except:
    pass

try:
    with open(pjoin(sys.prefix, 'hashstack_version.txt')) as f:
        hashstack = f.readline().strip()
except:
    pass

try:
    with open(pjoin(sys.prefix, 'proteus_version.txt')) as f:
        proteus = f.readline().strip()
    if os.path.exists(proteus):
        proteus = subprocess.check_output(git_cmd, cwd=proteus).strip()
except:
    pass
