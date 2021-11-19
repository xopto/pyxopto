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

ARG BASE_CONTAINER=pyxopto/jupyter

FROM ${BASE_CONTAINER}

LABEL maintainer="info@xopto.eu"

ARG TENSORFLOW_VERSION="tensorflow-gpu==2.3.*"
ARG PYTORCH_VERSION="torch==1.7.* torchvision"

USER root

# install gpu support
RUN python3 -m pip install --no-cache-dir pycuda

# install deep-learning
RUN python3 -m pip install --no-cache-dir $TENSORFLOW_VERSION
RUN python3 -m pip install --no-cache-dir $PYTORCH_VERSION

# switch to user
USER $NB_USER
