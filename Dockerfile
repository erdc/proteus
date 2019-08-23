FROM erdc/stack_base:latest

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
    make N=4 develop && \
    CC=gcc CXX=g++ ./linux/bin/pip3 install matplotlib && \
    PATH=/home/$NB_USER/proteus/linux/bin:$PATH make jupyter && \
    ./linux/bin/pip3 install jupyterhub && \
    rm -rf build && \
    rm -rf air-water-vv && \
    rm -rf .git && \
    rm -rf stack/.git && \
    rm -rf /home/$NB_USER/.cache 

ENV PATH /home/$NB_USER/proteus/linux/bin:$PATH
ENV LD_LIBRARY_PATH /home/$NB_USER/proteus/linux/lib:$LD_LIBRARY_PATH

USER root

CMD ["start-notebook.sh"]

# Add local files as late as possible to avoid cache busting
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/start.sh /usr/local/bin/start.sh
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/start-notebook.sh /usr/local/bin/start-notebook.sh
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/start-singleuser.sh /usr/local/bin/start-singleuser.sh

RUN chmod a+rx /usr/local/bin/*

RUN ipython kernel install

USER $NB_USER

# Import matplotlib the first time to build the font cache.
#ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
#RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot"
