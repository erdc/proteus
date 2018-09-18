#!/usr/bin/env python
import os
import optparse
import sys
import proteus

parser = optparse.OptionParser()
parser.add_option("-l", "--log",
                  help="Store information about what the code is doing, 0=none, 10=everything",
                  action="store",
                  type="int",
                  dest="logLevel",
                  default=1)
parser.add_option("--clean",
                  help="Clean folder",
                  action="store_true",
                  dest="clean")
parser.add_option("-C",
                  help="Context options",
                  action="store",
                  type="string",
                  default="",
                  dest="context")
parser.add_option("-D", "--dataDir",
                  help="Data directory",
                  action="store",
                  type="string",
                  dest="dataDir")
parser.add_option("-n","--num_proc",
                  help="Number of processors to use",
                  type="int",
                  action="store",
                  default=1,
                  dest="num_proc")

(opts,args) = parser.parse_args()

################
# CLEAN FOLDER #
################
if opts.clean is not None:
    os.system("rm -r *face *.csv __pycache__ *mesh* *.poly *.pyc *.log *.edge *.ele *.neig *.node *.h5 *.xmf *~ *#* *.txt *smb *pos *dmg splitMesh")
    exit()

############
# DATA DIR #
############
dataDir = "" if opts.dataDir is None else "-D " + opts.dataDir 

##############
# CALL PARUN #
##############
os.system("parun --SWEs " +
          "-l" + str(opts.logLevel) +
          " -v SWEs_so.py " +
          dataDir +
          " -C '" + opts.context + "'")
