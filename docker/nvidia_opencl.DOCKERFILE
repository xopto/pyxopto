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

ARG UBUNTU_VERSION=18.04
ARG CUDA_VERSION=11.0
ARG CUDA_TAG=devel
ARG CUDNN_VERSION=8

FROM nvidia/cuda:${CUDA_VERSION}-cudnn${CUDNN_VERSION}-${CUDA_TAG}-ubuntu${UBUNTU_VERSION}

LABEL maintainer="info@xopto.eu"

# switch to noninteractive apt-get mode
ARG DEBIAN_FRONTEND=noninteractive

# to avoid warnings about apt-utils
#RUN apt-get update && apt-get install --assume-yes apt-utils

# install opencl support
RUN apt-get update && apt-get install -y --no-install-recommends \
	ocl-icd-libopencl1 \
	clinfo \
	ocl-icd-opencl-dev

RUN rm -rf /var/lib/apt/lists/*

RUN mkdir -p /etc/OpenCL/vendors && \
	echo "libnvidia-opencl.so.1" > /etc/OpenCL/vendors/nvidia.icd

# nvidia-container-runtime
ENV NVIDIA_VISIBLE_DEVICES all
ENV NVIDIA_DRIVER_CAPABILITIES compute,utility
