#!/bin/sh

#EXE_PATH=/users/osahni/develop_bladapt/bladapt/blMesher
#EXE_PATH=/lore/zhanga/blMesher/
EXE_PATH=/lore/zhanga/blMesher/bladapt/blMesher/

export LD_LIBRARY_PATH=/net/common/meshSim/latest/lib/x64_rhel5_gcc41/psKrnl:$LD_LIBRARY_PATH

logfile=blMesher.log
#logfile=/dev/stdout

echo
echo "Running BL mesher (make sure to check $logfile)..."

#/usr/local/toolworks/totalview.8.11.0-0/bin/totalview /usr/local/openmpi/latest/bin/mpirun -a -np 1 $EXE_PATH/bin/x86_64_linux/BLMesher geom.x_t geom.sms 2 1 BLattr.inp -mesh_specify tetcoords1_2.dat > $logfile 2>&1

$EXE_PATH/bin/x86_64_linux/BLMesher Couette.x_t Splashcube.sms 2 1 BLattr.inp  > $logfile 2>&1
