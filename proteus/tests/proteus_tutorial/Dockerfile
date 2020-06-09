FROM erdc/proteus:gcc8
MAINTAINER Proteus Project <proteus@googlegroups.com>
USER root
COPY . /home/$NB_USER/proteus_tutorial
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
WORKDIR /home/$NB_USER/proteus_tutorial
RUN pip install --no-cache-dir notebook==5.*
