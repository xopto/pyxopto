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

. ./build_nvidia_jupyter.sh

docker build \
    --build-arg BASE_CONTAINER=xopto/pyxopto-nvidia-jupyter:$VERSION \
    --build-arg TENSORFLOW_VERSION='tensorflow-gpu==2.4.*' \
    --build-arg PYTORCH_VERSION='torch==1.7.1+cu110 torchvision==0.8.2+cu110 torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html' \
    -t xopto/pyxopto-nvidia-jupyter-dl:$VERSION \
    --file pyxopto_jupyter_dl.DOCKERFILE \
    ..


# Tensorflow CUDA compatibility matrix: https://www.tensorflow.org/install/source#gpu
#
# Run as 
#    "docker run --runtime=nvidia --rm -p 8888:8888 -it pyxopto/nvidia-jupyter-dl"
#
# Test tensorflow version < 1.12 in shell as:
#	python3 -c "import tensorflow as tf; tf.enable_eager_execution(); print(tf.reduce_sum(tf.random_normal([1000, 1000])))"
# Test tensorflow version > 1.12 in shell as:
#	python3 -c "import tensorflow as tf; print(tf.reduce_sum(tf.random.normal([1000, 1000])))"
