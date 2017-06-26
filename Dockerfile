FROM erdc/stack_base:latest

MAINTAINER Proteus Project <proteus@googlegroups.com>


USER jovyan

WORKDIR /home/$NB_USER

ENV CC mpicc
ENV CXX mpicxx
ENV F77 mpif77
ENV F90 mpif90

RUN cd proteus && make develop

ENV PATH /home/$NB_USER/proteus/linux2/bin:$PATH
ENV LD_LIBRARY_PATH /home/$NB_USER/proteus/linux2/lib:$LD_LIBRARY_PATH

RUN cd proteus && make jupyter

USER root

RUN jupyter kernelspec install-self

USER $NB_USER

# Import matplotlib the first time to build the font cache.
ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot"
