# Dockerfile for building preCICE on ubuntu 18.04

# Using ubuntu 18.04 as basis
FROM ubuntu:18.04

# Installing necessary dependacies for preCICE, boost 1.65 from apt-get
RUN apt-get -qq update && apt-get -qq install \
    build-essential \
    libboost-all-dev \
    libeigen3-dev \
    libxml2-dev \
    petsc-dev \
    git \
    python-numpy \
    python-dev \
    wget \
    bzip2 \
    cmake && \
    rm -rf /var/lib/apt/lists/*

# Rebuild image if force_rebuild after that command
ARG CACHEBUST

# create user precice
ARG uid=1000
ARG gid=1000
RUN groupadd -g ${gid} precice \
    && useradd -u ${uid} -g ${gid} -m -s /bin/bash precice
USER precice

# Setting some environment variables for installing preCICE
ENV CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:/usr/include/eigen3" \
    CPATH="/usr/include/eigen3:${CPATH}" \
    PETSC_DIR="/usr/lib/petscdir/3.6.2/" \
    PETSC_ARCH="x86_64-linux-gnu-real"

COPY . /home/precice/precice_develop

WORKDIR /home/precice/precice

RUN mkdir /home/precice/precice_build && \
    mkdir /home/precice/precice_install && \
    cd /home/precice/precice_build && \
    cmake -DBUILD_SHARED_LIBS=ON \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_INSTALL_PREFIX=/home/precice/precice-install \
          -DPRECICE_PETScMapping=yes \
          -DPRECICE_MPICommunication=yes \
          -DPRECICE_PythonActions=no \
          /home/precice/precice_develop && \
    make -j 10 && \
    make install && \
    ctest -V -R precice.parallel
  