#!/bin/bash -e
if [ $# -ne 3 ]; then 
	echo "Usage: $0 <first step> <last step> <num procs>"
else 
  first=$1
  last=$2
  procs=$3
  for step in `seq $first $last`; do
	set +e  # OK if ls fails
	x=`ls channel_p${step}0.h5`; 
	set -e
	if [ ! -z "$x" ]; then
      echo "step $step"
	   /users/behaff/proteusScorec/linux2-scorec/bin/gatherArchives.py -s $procs -f channel_p$step 
	else 
	  echo "skipping step $step"
	fi
  done
fi
