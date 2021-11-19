################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################

ARG BASE_CONTAINER=cuda-opencl-11.0

FROM ${BASE_CONTAINER}

LABEL maintainer="info@xopto.eu"

ARG NB_USER="jovyan"
ARG NB_UID="1000"
ARG NB_GID="100"

ARG TENSORFLOW_VERSION="tensorflow-gpu==2.3.*"
ARG PYTORCH_VERSION="torch==1.7.* torchvision"

ENV DEBIAN_FRONTEND noninteractive

# expose the jupyter port
EXPOSE 8888

RUN apt-get update && apt-get -yq dist-upgrade \
    tini \
    sudo \
    unzip \
    locales \
    git \
    texlive-latex-recommended \
    gnupg2 \
    nano \
    gmsh \
    gcc \
    nodejs \
    libboost-all-dev \
    python3 \
    python3-pip

# Configure environment
ENV SHELL=/bin/bash \
    NB_USER=$NB_USER \
    NB_UID=$NB_UID \
    NB_GID=$NB_GID \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    HOME=/home/$NB_USER \
    PYTHONPATH=/home/$NB_USER/pyxopto \
    PYXOPTO_USER_PATH=/home/$NB_USER/work/.xopto/pyxopto

# install python packages
RUN python3 -m pip install --no-cache-dir --upgrade pip
RUN python3 -m pip install --no-cache-dir setuptools
RUN python3 -m pip install --no-cache-dir wheel
RUN python3 -m pip install --no-cache-dir numpy
RUN python3 -m pip install --no-cache-dir scipy
RUN python3 -m pip install --no-cache-dir matplotlib
RUN python3 -m pip install --no-cache-dir jinja2
RUN python3 -m pip install --no-cache-dir shapely
RUN python3 -m pip install --no-cache-dir numba
RUN python3 -m pip install --no-cache-dir pybind11
RUN python3 -m pip install --no-cache-dir python-scattnlay
RUN python3 -m pip install --no-cache-dir pandas
RUN python3 -m pip install --no-cache-dir keras

# install gpu support
RUN python3 -m pip install --no-cache-dir pyopencl

# install jupyter
RUN python3 -m pip install --no-cache-dir jupyter
RUN python3 -m pip install --no-cache-dir jupyterlab

COPY docker/script/start_jupyter.sh /start_jupyter.sh

# Configure container entrypoint
ENTRYPOINT ["tini", "-g", "--"]

# Start jupyter server by default
CMD ["/bin/bash", "-c", "/start_jupyter.sh"]

RUN echo "LC_ALL=en_US.UTF-8" >> /etc/environment
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen
RUN echo "LANG=en_US.UTF-8" > /etc/locale.conf
RUN locale-gen en_US.UTF-8

# create a new user ... delete home ... just in case if exists
RUN rm -rf /home/$NB_USER
RUN useradd -m -s /bin/bash -G video,sudo -N -u $NB_UID $NB_USER
RUN chown $NB_USER:$NB_GID /home/$NB_USER && \
	chmod 755 /home/$NB_USER

# switch to user
USER $NB_USER
WORKDIR $HOME

# create the user work directory
RUN mkdir -p ~/work

# create the user pyxopto directory
RUN mkdir -p ~/pyxopto/xopto
RUN mkdir -p ~/pyxopto/examples
RUN mkdir -p ~/pyxopto/test
RUN mkdir -p ~/pyxopto/documentation

# copy pyxopto into the container
COPY xopto                   $HOME/pyxopto/xopto/
COPY examples                $HOME/pyxopto/examples/
COPY test                    $HOME/pyxopto/test/
COPY docs/build/docs/html    $HOME/pyxopto/documentation/
