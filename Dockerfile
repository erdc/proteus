FROM erdc/stack_base:gcc8

MAINTAINER Proteus Project <proteus@googlegroups.com>

USER root

USER jovyan

WORKDIR /home/$NB_USER

RUN rm -rf proteus && \
    git clone https://github.com/erdc/proteus && \
    cd proteus && \
    git checkout check_docker_hub && \
    git submodule update --init --recursive && \
    gcc --version && \
    g++ --version && \
    gfortran --version && \
    ulimit -m 2000000 && \
    ulimit -v 4000000 && \
    ulimit -a && \
    free -m && \
    make PROTEUS_OPT="-DNDEBUG -g0 -O0" N=1 develop && \
    rm -rf air-water-vv && \
    rm -rf .git && \
    rm -rf stack/.git && \
    rm -rf /home/$NB_USER/.cache 

ENV PATH /home/$NB_USER/proteus/linux/bin:$PATH
ENV LD_LIBRARY_PATH /home/$NB_USER/proteus/linux/lib:/home/$NB_USER/proteus/linux/lib64:$LD_LIBRARY_PATH
