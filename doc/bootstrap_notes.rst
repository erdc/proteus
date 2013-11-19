Installing Proteus and its dependencies on a SuSE 11.1 VM
=========================================================

Downloading SuSE
----------------

You'll need to log in to a Novell account, but the file itself is free.

You are looking for ``SLES-11-SP1-DVD-x86_64-GM-DVD1.iso`` with a
checksum of ``d2e10420f3689faa49a004b60fb396b7``

SuSE Setup
----------

For the below requirements, these instructions install the majority of packages into
``/usr/local``. If you do not have access to an account with
administrative privileges, you could try a local build, but note that it
is going to be very difficult for you to install the GFortran compiler.

Manual - Firewall/NFSServer/SSHD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

YaST and zypper are your friends.

Fortran
~~~~~~~

sudo zypper addrepo
http://download.opensuse.org/repositories/home:zhy20120210:SLES-11-SP1-x86-64/SLE\_11/home:zhy20120210:SLES-11-SP1-x86-64.repo
sudo zypper install gcc43-fortran sudo zypper install gcc-fortran

**NOTE**: This renders your system as vendor-unsupported because the
official SLES distribution does not contain a ``gfortran`` binary. The
major/minor GNU Compiler toolchain numbers are matched, but this will
give you a slightly different toolchain than what SLES officially
distributes and supports.

OpenSSL
~~~~~~~

A variety of current tools require a modern OpenSSL.

::

    wget http://www.openssl.org/source/openssl-1.0.1e.tar.gz
    tar -zxvf openssl-1.0.1e.tar.gz
    cd openssl-1.0.1e
    ./config shared zlib-dynamic
    make
    sudo make install
    sudo mv /usr/local/lib64/libcrypto.a{,.backup}
    sudo mv /usr/local/lib64/libssl.a{,.backup}
    sudo ln -s /usr/local/ssl/lib/libssl.so /lib64
    sudo ln -s /usr/local/ssl/lib/libcrypto.so /lib64
    sudo ln -s /usr/local/ssl/lib/libssl.so.1.0.0 /lib64
    sudo ln -s /usr/local/ssl/lib/libcrypto.so.1.0.0 /lib64
    cd ..

Curl
~~~~

::

    wget http://curl.haxx.se/download/curl-7.33.0.tar.gz
    tar -zxvf curl-7.33.0.tar.gz

    cd curl-7.33.0
    ./configure --with-ssl --prefix=/usr/local LDFLAGS="-ldl"
    make
    sudo make install
    cd ..

Expat
~~~~~

::

    wget http://sourceforge.net/projects/expat/files/latest/download
    tar -zxvf expat-2.1.0.tar.gz
    cd expat-2.1.0
    ./configure --prefix=/usr/local
    make
    sudo make install
    cd ..

Git
~~~

::

    wget https://github.com/git/git/archive/v1.8.4.1.tar.gz
    tar -zxvf git.tar.gz
    cd git-1.8.4.1
    make
    sudo make install

MPICH
~~~~~

::

    sudo zypper install -y mpich2
    wget http://www.mpich.org/static/downloads/1.5/mpich2-1.5.tar.gz
    tar -zxvf mpich2-1.5.tar.gz
    cd mpich2-1.5
    ./configure --prefix=/usr/local --with-shared
    make
    sudo make install

You may also need to edit ``/etc/hosts`` and add an alias for 127.0.01
equal to the output of ``hostname``

Proteus
-------

::

    mkdir ~/old
    cd ~/old
    git clone git@github.com:erdc-cm/proteus.git
    cd proteus
    export PROTEUS=$(pwd)
    source envConfig/linux-suse.bash
    make

You should be able to test the build by heading to
``proteusModule/test`` and running
``test_meshPartitionFromTetgenFiles.py`` from there.
