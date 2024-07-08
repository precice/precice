# Dockerfile to build a ubuntu image containing the installed develop version of preCICE

FROM ubuntu:22.04
# Add the precice user
RUN useradd -m -s /bin/bash precice
ENV TZ=Europe/Berlin
RUN apt-get update && \
    apt-get -yy upgrade && \
    echo $TZ > /etc/timezone && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    apt-get -yy install git build-essential lsb-release cmake libeigen3-dev libxml2-dev libboost-all-dev python3-dev python3-numpy petsc-dev && \
    rm -rf /var/lib/apt/lists/*
# We do not use --depth=1 as we need the history for git describe in the version
RUN git clone --branch=develop https://github.com/precice/precice.git /root/precice && \
    cd /root/precice && \
    cmake --preset=debian-package && \
    cd build && \
    make -j $(nproc) && \
    cpack -G DEB && \
    apt-get -yy install ./libprecice*.deb && \
    rm -rf /root/precice
# Make sure the installation is functional
RUN precice-tools version
