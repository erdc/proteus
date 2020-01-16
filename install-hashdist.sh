#!/bin/bash

git lfs pull
./stack/hit/bin/hit init-home
./stack/hit/bin/hit remote add http://levant.hrwallingford.com/hashdist_src --objects="source"
./stack/hit/bin/hit remote add http://levant.hrwallingford.com/hashdist_ubuntu_16_04 --objects="build"
make stack/default.yaml
pushd $HOME
mkdir -p hashdist_src
mkdir -p hashdist_bld
rm -rf .hashdist/src .hashdist/bld
ln -s $HOME/hashdist_src .hashdist/src
ln -s $HOME/hashdist_bld .hashdist/bld
popd
pushd stack
echo $PWD
ls -l
./hit/bin/hit build -j 2 -v default.yaml
popd
export PATHSAVE=$PATH
export PATH=$PWD/linux/bin:$PATH
export LD_LIBRARY_PATH=$PWD/linux/lib:$LD_LIBRARY_PATH
PROTEUS_OPT="-w -O2 -UNDEBUG" FC=gfortran CC=mpicc CXX=mpicxx make develop N=2
export SSL_CERT_DIR=/etc/ssl/certs
#./linux/bin/pip3 install matplotlib
