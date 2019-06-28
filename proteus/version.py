"""
Determine versions of proteus and stack by inspecting git repository
"""

import os
from os.path import join as pjoin
import sys
import subprocess

proteus_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            os.pardir)
stack_dir = os.path.join(proteus_dir,'stack')

git_cmd = ['git', 'log', '-1', '--pretty=%H']

stack = 'unknown'
proteus = 'unknown'

try:
    stack = str(subprocess.check_output(git_cmd, cwd=stack_dir).strip(),'utf-8')
except:
    pass

try:
    proteus = str(subprocess.check_output(git_cmd, cwd=proteus_dir).strip(),'utf-8')
except:
    pass
