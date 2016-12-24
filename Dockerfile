FROM debian:jessie

MAINTAINER Proteus Project <proteus@googlegroups.com>

USER root

# Install all OS dependencies for fully functional notebook server
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -yq --no-install-recommends --fix-missing \
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
    freeglut3-dev \
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
    && apt-get clean

RUN pip3 install notebook terminado

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Install Tini
RUN wget --quiet https://github.com/krallin/tini/releases/download/v0.9.0/tini && \
    echo "faafbfb5b079303691a939a747d7f60591f2143164093727e870b289a44d9872 *tini" | sha256sum -c - && \
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

RUN mkdir /home/$NB_USER/work && \
    mkdir /home/$NB_USER/.jupyter && \
    mkdir /home/$NB_USER/.local && \
    mkdir /home/$NB_USER/.hashdist 

ADD https://dl.dropboxusercontent.com/u/26353144/hashdist_config_jovyan.yaml /home/$NB_USER/.hashdist/config.yaml

RUN chown -R $NB_USER:users /home/$NB_USER

USER jovyan

RUN ls /home/$NB_USER/.hashdist && \
    cat /home/$NB_USER/.hashdist/config.yaml

WORKDIR /home/$NB_USER

RUN git clone https://github.com/erdc-cm/workshops -b erdc-fsi-tutorials

RUN cat /home/$NB_USER/.hashdist/config.yaml && \
    git clone https://github.com/erdc-cm/proteus && \
    cd proteus && git checkout update_jupyter && \
    make hashdist stack stack/default.yaml && \
    cd stack && \
    /usr/bin/mpicc -show && \
    which mpiexec && \
    which mpirun && \
    ../hashdist/bin/hit build default.yaml -v

ENV CC mpicc
ENV CXX mpicxx
ENV F77 mpif77
ENV F90 mpif90

RUN cd proteus && make develop

ENV PATH /home/$NB_USER/proteus/linux2/bin:$PATH
ENV LD_LIBRARY_PATH /home/$NB_USER/proteus/linux2/lib:$LD_LIBRARY_PATH

RUN make jupyter

USER root

# Configure container startup as root
EXPOSE 8888
WORKDIR /home/$NB_USER/work
ENTRYPOINT ["tini", "--"]
CMD ["start-notebook.sh"]

RUN cd /usr/local/bin && \
    wget https://raw.githubusercontent.com/jupyter/docker-stacks/master/minimal-notebook/start-notebook.sh

ADD https://raw.githubusercontent.com/jupyter/docker-stacks/master/minimal-notebook/jupyter_notebook_config.py /home/$NB_USER/.jupyter/jupyter_notebook_config.py

RUN mkdir /etc/jupyter && \
    chmod a+rwX /etc/jupyter && \
    chown -R $NB_USER:users /home/$NB_USER

#jupyter/ipython extensions
RUN pip3 install \
    ipyparallel \
    ipywidgets

# Switch back to jovyan to avoid accidental container runs as root
USER jovyan
    
RUN cd ~/.jupyter && \
    ipython profile create mpi --parallel && \
    ipcluster nbextension enable && \
    echo '\nc.NotebookApp.server_extensions.append("ipyparallel.nbextension")' >> /home/$NB_USER/.jupyter/jupyter_notebook_config.py && \
    cp jupyter_notebook_config.py /etc/jupyter/ && \
    echo "c.LocalControllerLauncher.controller_cmd = ['python2', '-m', 'ipyparallel.controller']\nc.LocalEngineSetLauncher.engine_cmd = ['python2', '-m', 'ipyparallel.engine']\n" \
          >> /home/$NB_USER/.ipython/profile_mpi/ipcluster_config.py  

USER root

RUN jupyter kernelspec install-self

# fetch juptyerhub-singleuser entrypoint
RUN wget -q https://raw.githubusercontent.com/jupyter/jupyterhub/master/scripts/jupyterhub-singleuser -O /usr/local/bin/jupyterhub-singleuser && \
    chmod 755 /usr/local/bin/jupyterhub-singleuser

ADD https://raw.githubusercontent.com/jupyter/dockerspawner/master/singleuser/singleuser.sh /srv/singleuser/singleuser.sh

RUN chmod 755 /srv/singleuser/singleuser.sh

USER jovyan

RUN cat /srv/singleuser/singleuser.sh
# smoke test that it's importable at least
RUN sh /srv/singleuser/singleuser.sh -h
CMD ["sh", "/srv/singleuser/singleuser.sh"]
