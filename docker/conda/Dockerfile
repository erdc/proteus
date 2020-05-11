FROM erdc/conda_base:latest

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
    pip install -v -e . && \ 
    rm -rf build && \
    rm -rf air-water-vv && \
    rm -rf .git && \
    rm -rf stack/.git && \
    rm -rf /home/$NB_USER/.cache 

