
from subprocess import call
import os.path

source_dir = "/opt/proteus/"
target_dir = "/usr/local/work/"


if (not os.path.islink(target_dir + "stack")):
    call(["ln", "-s", source_dir + "stack",  target_dir + "stack"])

if (not os.path.islink(target_dir + "hashdist")):
    call(["ln", "-s", source_dir + "hashdist",  target_dir + "hashdist"])

if (not os.path.islink(target_dir + "linux2")):
    call(["ln", "-s", source_dir + "linux2",  target_dir + "linux2"])

