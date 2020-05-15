FROM erdc/stack_base:gcc9

MAINTAINER Proteus Project <proteus@googlegroups.com>

USER jovyan

WORKDIR /home/$NB_USER

ENV CC mpicc
ENV CXX mpicxx
ENV F77 mpif77
ENV F90 mpif90

RUN rm -rf proteus && \
    git clone https://github.com/erdc/proteus && \
    cd proteus && \
    git checkout master && \
    git submodule update --init --recursive && \
    gcc --version && \
    g++ --version && \
    gfortran --version && \
    make PROTEUS_OPT="-DNDEBUG -g0 -O0" N=1 develop && \
    . ./linux/bin/proteus_env.sh \
    mpicc --version && \
    mpicxx --version && \
    mpifort --version && \
    rm -rf air-water-vv && \
    rm -rf .git && \
    rm -rf stack/.git && \
    rm -rf /home/$NB_USER/.cache 

ENV PATH /home/$NB_USER/proteus/linux/bin:$PATH
ENV LD_LIBRARY_PATH /home/$NB_USER/proteus/linux/lib:/home/$NB_USER/proteus/linux/lib64:$LD_LIBRARY_PATH
