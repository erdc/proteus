#
#garent.pgi
#
module load acml
module load gcc #load gcc compiler in addition to pgi
export PROTEUS_ARCH=garnet.pgi
export PROTEUS_PREFIX=${PROTEUS}/${PROTEUS_ARCH}
export PROTEUS_PYTHON=${PROTEUS_PREFIX}/bin/python
export PATH=.:${PROTEUS_PREFIX}/bin:${HOME}/src/proteus/garnet/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PROTEUS_PREFIX}/lib:${ACML_DIR}/pgi64/lib:/opt/cray/udreg/2.2-1.0301.2966.16.2.gem/lib64:/opt/cray/ugni/2.1-1.0301.2967.10.23.gem/lib64:/opt/cray/pmi/1.0-1.0000.8160.39.2.gem/lib64:/opt/cray/dmapp/3.0-1.0301.2968.22.24.gem/lib64:/opt/cray/xpmem/0.1-2.0301.25333.20.2.gem/lib64:/opt/cray/mpt/5.1.2/xt/gemini/mpich2-pgi/lib:/opt/acml/4.4.0/pgi64/lib:/opt/xt-libsci/10.4.9/pgi/lib:/opt/cray/portals/default/lib64:/opt/cray/xe-sysroot/3.1.61.securitypatch.20110308/usr/lib64:/opt/cray/xe-sysroot/3.1.61.securitypatch.20110308/lib64:/opt/cray/xe-sysroot/3.1.61.securitypatch.20110308:/usr/lib/alps:/usr/lib/alps:/opt/pgi/10.9.0/linux86-64/10.9/libso:/opt/pgi/10.9.0/linux86-64/10.9/lib:/usr/lib64/gcc/x86_64-suse-linux/4.3
