#!/bin/bash

export MPLBACKEND="AGG"
PATH=./linux/bin:$PATH py.test -n 1 --dist=loadfile --forked -v proteus/tests --ignore proteus/tests/POD  --ignore proteus/tests/MeshAdaptPUMI --cov=proteus
