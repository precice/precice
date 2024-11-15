# Dockerfile to build a ubuntu image containing the installed Debian package of a release

FROM ubuntu:22.04
# Add the precice user
RUN useradd -m -s /bin/bash precice
# Fix the installation of tzdata for Ubuntu
ARG TIMEZONE=Europe/Berlin
RUN export TZ=$TIMEZONE && echo $TZ > /etc/timezone && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    apt-get -yy update && apt-get -yy install wget tzdata lsb-release && rm -rf /var/lib/apt/lists/*
# The version to fetch the package for: X.Y.Z
ARG version
RUN echo "$version" | grep "[0-9]\+\.[0-9]\+\.[0-9]\+" > /dev/null # The version must follow the format X.Y.Z
RUN wget -q -O libprecice.deb https://github.com/precice/precice/releases/download/v${version}/libprecice`echo ${version} | sed 's/\([0-9]\+\)\.\([0-9]\+\.[0-9]\+\)/\1_\1.\2/'`_$(lsb_release -sc).deb && \
    apt-get update && apt-get -yy install ./libprecice.deb && \
    rm libprecice.deb && rm -rf /var/lib/apt/lists/*
# Make sure the installation is functional
RUN precice-tools version
