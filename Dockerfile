FROM debian:jessie

MAINTAINER Proteus Project <proteus@googlegroups.com>

USER root

ENV DEBIAN_FRONTEND noninteractive

RUN REPO=http://cdn-fastly.deb.debian.org \
    && echo "deb $REPO/debian jessie main\ndeb $REPO/debian-security jessie/updates main" > /etc/apt/sources.list \
    && apt-get update && apt-get -yq dist-upgrade \
    && apt-get install -yq --no-install-recommends --fix-missing \
    git \
    vim \
    jed \
    emacs \
    wget \
    build-essential \
    python-dev \
    ca-certificates \
    bzip2 \
    unzip \
    libsm6 \
    pandoc \
    texlive-latex-base \
    texlive-latex-extra \
    texlive-fonts-extra \
    texlive-fonts-recommended \
    texlive-generic-recommended \
    sudo \
    locales \
    libxrender1 \
    libav-tools \
    libmpich2-dev \
    liblapack-dev \
    freeglut3 \
    freeglut3-dev \
    libglew1.5 \
    libglew1.5-dev \
    libglu1-mesa \
    libglu1-mesa-dev \
    libgl1-mesa-glx \
    libgl1-mesa-dev \
    curl \
    libjpeg-dev \
    m4 \
    libssl-dev \
    ssh \
    mpich2 \
    python3 \
    python3-pip \
    python3-doc \
    python3-tk \
    python3-venv \
    python3-genshi \
    python3-lxml \
    python3-openssl \
    python3-pyasn1 \
    python3.4-venv \
    python3.4-doc \
    binfmt-support \
    python3-dev \
    python3-wheel \
    libffi-dev \
    python-lzma \
    python-pip \
    cmake \
    gfortran \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Install Tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.10.0/tini && \
    echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - && \
    mv tini /usr/local/bin/tini && \
    chmod +x /usr/local/bin/tini

# Configure environment
ENV SHELL /bin/bash
ENV NB_USER jovyan
ENV NB_UID 1000
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Create jovyan user with UID=1000 and in the 'users' group
RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER

RUN mkdir /home/$NB_USER/.jupyter && \
    mkdir /home/$NB_USER/.local && \
    echo "cacert=/etc/ssl/certs/ca-certificates.crt" > /home/$NB_USER/.curlrc

RUN chown -R $NB_USER:users /home/$NB_USER

USER jovyan

WORKDIR /home/$NB_USER

RUN git clone https://github.com/erdc-cm/proteus && \
    cd proteus && \
    make hashdist stack stack/default.yaml && \
    ./hashdist/bin/hit init-home && \
    ./hashdist/bin/hit remote add https://dl.dropboxusercontent.com/u/26353144/hashdist_src --objects="source" && \
    ./hashdist/bin/hit remote add https://dl.dropboxusercontent.com/u/26353144/hashdist_jessie --objects="build" && \
    cd stack && \
    ../hashdist/bin/hit build default.yaml -v

ENV CC mpicc
ENV CXX mpicxx
ENV F77 mpif77
ENV F90 mpif90

RUN cd proteus && make develop

ENV PATH /home/$NB_USER/proteus/linux2/bin:$PATH
ENV LD_LIBRARY_PATH /home/$NB_USER/proteus/linux2/lib:$LD_LIBRARY_PATH

RUN cd proteus && make jupyter

USER root

RUN pip3 install pyzmq==16.0.2 --install-option="--zmq=/home/$NB_USER/proteus/linux2"
RUN pip3 install ipyparallel==6.0.2 ipython==5.3.0 terminado==0.6 jupyter==1.0.0 jupyterlab==0.18.1  ipywidgets==6.0.0 ipyleaflet==0.3.0 jupyter_dashboards==0.7.0 pythreejs==0.3.0 rise==4.0.0b1 cesiumpy==0.3.3 bqplot==0.9.0
RUN /usr/local/bin/jupyter serverextension enable --py jupyterlab --sys-prefix \
    && /usr/local/bin/jupyter nbextension enable --py --sys-prefix widgetsnbextension \
    && /usr/local/bin/jupyter nbextension enable --py --sys-prefix bqplot \
    && /usr/local/bin/jupyter nbextension enable --py --sys-prefix pythreejs \
    && /usr/local/bin/jupyter nbextension enable --py --sys-prefix ipyleaflet \
    && /usr/local/bin/jupyter nbextension install --py --sys-prefix rise \
    && /usr/local/bin/jupyter nbextension enable --py --sys-prefix rise \
    && /usr/local/bin/jupyter dashboards quick-setup --sys-prefix \
    && /usr/local/bin/jupyter nbextension install --sys-prefix --py ipyparallel \
    && /usr/local/bin/jupyter nbextension enable --sys-prefix --py ipyparallel \
    && /usr/local/bin/jupyter serverextension enable --sys-prefix --py ipyparallel

EXPOSE 8888
WORKDIR /home/$NB_USER

ENTRYPOINT ["tini", "--"]
CMD ["start-notebook.sh"]

# Add local files as late as possible to avoid cache busting
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/start.sh /usr/local/bin/start.sh
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/start-notebook.sh /usr/local/bin/start-notebook.sh
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/start-singleuser.sh /usr/local/bin/start-singleuser.sh
ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/base-notebook/jupyter_notebook_config.py /home/$NB_USER/.jupyter/jupyter_notebook_config.py

RUN chmod a+rx /usr/local/bin/*

RUN chown -R $NB_USER:users /home/$NB_USER/.jupyter

RUN jupyter kernelspec install-self

# Switch back to jovyan to avoid accidental container runs as root
USER $NB_USER

RUN cd ~/.jupyter && \
    ipython profile create mpi --parallel && \
    ipcluster nbextension enable --user && \
    echo '\nc.NotebookApp.server_extensions.append("ipyparallel.nbextension")' >> /home/$NB_USER/.jupyter/jupyter_notebook_config.py && \
    echo "c.LocalControllerLauncher.controller_cmd = ['python2', '-m', 'ipyparallel.controller']\nc.LocalEngineSetLauncher.engine_cmd = ['python2', '-m', 'ipyparallel.engine']\n" \
          >> /home/$NB_USER/.ipython/profile_mpi/ipcluster_config.py \
    && jupyter serverextension enable --py jupyterlab --user \
    && jupyter nbextension enable --py --user widgetsnbextension\
    && jupyter nbextension enable --py --user bqplot \
    && jupyter nbextension enable --py --user pythreejs \
    && jupyter nbextension enable --py --user ipyleaflet \
    && jupyter nbextension install --py --user rise \
    && jupyter nbextension enable --py --user rise

# Import matplotlib the first time to build the font cache.
ENV XDG_CACHE_HOME /home/$NB_USER/.cache/
RUN MPLBACKEND=Agg python -c "import matplotlib.pyplot"
