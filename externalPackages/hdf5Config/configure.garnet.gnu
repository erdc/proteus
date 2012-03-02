LDFLAGS="-L${PROTEUS_PREFIX}/lib" LIBS="-lsz -lz" CFLAGS="-I${PROTEUS_PREFIX}/include" ./configure \
--with-szlib=${PROTEUS_PREFIX} --with-zlib=${PROTEUS_PREFIX} \
--enable-shared --enable-threadsafe -disable-parallel --disable-fortran --disable-cxx --with-pthread --with-pic \
--prefix=${PROTEUS_PREFIX}
