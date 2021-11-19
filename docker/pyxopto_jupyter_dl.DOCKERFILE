ARG BASE_CONTAINER=pyxopto/jupyter

FROM ${BASE_CONTAINER}

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
