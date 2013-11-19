x=$(patch -Nup1 < ../mpi4pyConfig/no-python-mpi.patch.garnet)
echo "$x"
${PROTEUS_PYTHON} setup.py build build_exe --mpicxx=CC --mpicc=cc --mpif77=ftn --mpif90=ftn --configure
