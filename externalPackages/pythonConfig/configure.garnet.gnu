export LD_LIBRARY_PATH=/usr/local/usp/openssl/lib:$LD_LIBRARY_PATH
./configure --prefix=${PROTEUS_PREFIX} --enable-shared
x=$(patch -Nup2 < ../pythonConfig/add_ssl_to_setup.patch.garnet)
echo $x
x=$(patch -Nup2 < ../pythonConfig/add_usr_local_usp_openssl_to_setup_py.patch.garnet)
echo $x
