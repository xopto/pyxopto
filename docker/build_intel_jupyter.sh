#!/bin/bash
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

set -e

export PYTHONPATH=$PWD/..

if [ -z "$VERSION" ]; then
    VERSION=v$(python3 -c "import xopto; print(xopto.__version__)")
fi

docker build \
    --build-arg BASE_CONTAINER=intelopencl/intel-opencl:ubuntu-20.04-ppa \
    -t xopto/pyxopto-intel-jupyter-base \
    --file pyxopto_jupyter.DOCKERFILE \
    ..

docker build \
    --build-arg BASE_CONTAINER=xopto/pyxopto-intel-jupyter-base \
    -t xopto/pyxopto-intel-jupyter:$VERSION \
    --file pyxopto_jupyter_finalize.DOCKERFILE \
    ..

# Tensorflow CUDA compatibility matrix: https://www.tensorflow.org/install/source#gpu
#
# Run as:
#   "docker run --rm -p 8888:8888 -it xopto/pyxopto-intel-jupyter:$VERSION"
