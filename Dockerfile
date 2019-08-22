FROM erdc/stack_base:latest

MAINTAINER Proteus Project <proteus@googlegroups.com>

USER jovyan

WORKDIR /home/$NB_USER

ENV CC mpicc
ENV CXX mpicxx
ENV F77 mpif77
ENV F90 mpif90

RUN cd proteus && git checkout master && git pull && make N=4 develop
RUN cd proteus && CC=gcc CXX=g++ ./linux/bin/pip install matplotlib

ENV PATH /home/$NB_USER/proteus/linux/bin:$PATH

RUN cd proteus && PATH=/usr/bin:/usr/local/bin:$PATH make jupyter

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
