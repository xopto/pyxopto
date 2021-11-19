ARG UBUNTU_VERSION=18.04
ARG CUDA_VERSION=11.0
ARG CUDA_TAG=devel
ARG CUDNN_VERSION=8

FROM nvidia/cuda:${CUDA_VERSION}-cudnn${CUDNN_VERSION}-${CUDA_TAG}-ubuntu${UBUNTU_VERSION}

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
