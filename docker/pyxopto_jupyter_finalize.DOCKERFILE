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

ARG BASE_CONTAINER

FROM ${BASE_CONTAINER}

LABEL maintainer="info@xopto.eu"

# switch to user
USER $NB_USER

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
