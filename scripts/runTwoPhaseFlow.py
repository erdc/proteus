#!/usr/bin/env python
import os
import optparse
import sys
import proteus

parser = optparse.OptionParser()
parser.add_option("-f", "--file_name",
                  help="Name of file to run",
                  action="store",
                  dest="file_name",
                  default=None)
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
parser.add_option("-n","--num_proc",
                  help="Number of processors to use",
                  type="int",
                  action="store",
                  default=1,
                  dest="num_proc")
parser.add_option("-C",
                  help="Context options",
                  action="store",
                  dest="context")
parser.add_option("-D", "--dataDir",
                  help="Data directory",
                  action="store",
                  type="string",
                  dest="dataDir")

# NAME OF PROBLEM #
#assert opts.file_name is not None, "Provide name of file to run: -f name_of_file.py"
#current_path = os.getcwd()

# PUMI
(opts,args) = parser.parse_args()

# SOME ASSERTS #
if opts.context is not None:
    assert "useSuperlu" not in opts.context, "Don't assign value to useSuperlu in -C"
    assert "usePUMI" not in opts.context, "Don't assign value to usePUMI in -C"

################
# CLEAN FOLDER #
################
if opts.clean is not None:
    os.system("rm -r *face *.csv __pycache__ *mesh* *.poly *.pyc *.log *.edge *.ele *.neig *.node *.h5 *.xmf *~ *#* *.txt *smb")
    exit()

#################
# PARALLEL RUNS #
#################
assert opts.num_proc >= 1, "Argument of -n must be an integer larger or equal to 1"
if opts.num_proc==1:
    useSuperlu="True"
else:
    useSuperlu="False"

###########
# CONTEXT #
###########
context = "" if opts.context is None else " "+opts.context

############
# DATA DIR #
############
dataDir = "" if opts.dataDir is None else "-D " + opts.dataDir 

##############
# CALL PARUN #
##############
os.system("mpirun -np " + str(opts.num_proc) +
          " parun --TwoPhaseFlow -l" + str(opts.logLevel) +
          " -v TwoPhaseFlow_so.py " +
          dataDir +
          " -C 'useSuperlu=" + useSuperlu + context + "'") 
