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
    --build-arg CUDA_VERSION=11.0.3 \
    --build-arg CUDNN_VERSION=8 \
    --build-arg UBUNTU_VERSION=20.04 \
    -t xopto/pyxopto-opencl-cuda-11.0.3 \
    --file nvidia_opencl.DOCKERFILE \
    .

docker build \
    --build-arg BASE_CONTAINER=xopto/pyxopto-opencl-cuda-11.0.3 \
	-t xopto/pyxopto-nvidia-jupyter:$VERSION \
    --file pyxopto_jupyter.DOCKERFILE \
    ..

# Tensorflow CUDA compatibility matrix: https://www.tensorflow.org/install/source#gpu
#
# Run as 
#    "docker run --runtime=nvidia --rm -p 8888:8888 -it pyxopto/nvidia-jupyter"
#
# Test tensorflow version < 1.12 in shell as:
#	python3 -c "import tensorflow as tf; tf.enable_eager_execution(); print(tf.reduce_sum(tf.random_normal([1000, 1000])))"
# Test tensorflow version > 1.12 in shell as:
#	python3 -c "import tensorflow as tf; print(tf.reduce_sum(tf.random.normal([1000, 1000])))"
