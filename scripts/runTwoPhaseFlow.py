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
parser.add_option("--pumi",
                  help="Run with adaptivity",
                  action="store_true",
                  default=False,
                  dest="pumi")
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
parser.add_option("-f","--fileName",
                  help="Name of setup file",
                  action="store",
                  type="string",
                  dest="fileName",
                  default="")

(opts,args) = parser.parse_args()

################
# CLEAN FOLDER #
################
if opts.clean is not None:
    os.system("rm -r *face *.csv __pycache__ *mesh* *.poly *.pyc *.log *.edge *.ele *.neig *.node *.h5 *.xmf *~ *#* *.txt *smb *pos *dmg splitMesh *pumi*")
    exit()

############
# DATA DIR #
############
dataDir = "" if opts.dataDir is None else "-D " + opts.dataDir

########
# PUMI #
########
if opts.pumi:
    from os import listdir
    from os.path import isfile, join
    path = "./splitMesh/"

    # FIRST STAGE: CREATE THE MESH IN SERIAL #

    err = os.system("parun --TwoPhaseFlow --pumi -l" + str(opts.logLevel) + " -v TwoPhaseFlow_so.py ")
    if(err):
        import sys
        sys.exit("error during proteus mesh to PUMI model derivation phase\n")

    # SECOND STAGE: SPLIT THE MESH IN PARALLEL #
    split = ("mpiexec -np " + str(opts.num_proc) +
             " split Reconstructed.dmg Reconstructed.smb " + path + " " + str(opts.num_proc))
    os.system("rm -r " + path)
    os.system(split)
    # rename files
    files = [f for f in listdir(path) if isfile(join(path, f))]
    [os.rename(path+f,path+"splitMesh"+str(counter)+".smb") for counter,f in enumerate(files)]

#############
# FILE NAME #
#############
fileName = ""
if opts.fileName != "":
    fileName = " -f " + opts.fileName

##############
# CALL PARUN #
##############
assert opts.num_proc >= 1, "Argument of -n must be an integer larger or equal to 1"
usePumi = "--pumi --pumiStage 2 " if opts.pumi else " "
if opts.num_proc==1:
    os.system("parun --TwoPhaseFlow " +
              usePumi +
              "-l" + str(opts.logLevel) +
              " -v TwoPhaseFlow_so.py " +
              dataDir +
              " -C '" + opts.context + "'" +
              fileName)
else:
    path_utils = proteus.__path__[0]+"/TwoPhaseFlow/utils/"
    petsc_options = path_utils + "petsc.options.asm"
    os.system("mpirun -np " + str(opts.num_proc) +
              " parun --TwoPhaseFlow --TpFlowParallel " +
              usePumi +
              "-l" + str(opts.logLevel) +
              " -v TwoPhaseFlow_so.py " +
              dataDir +
              " -O " + petsc_options +
              " -C '" + opts.context + "'" +
              fileName)

