# Dockerfile to build a ubuntu image containing the installed Debian package of a release 

FROM ubuntu:20.04
# Fix the installation of tzdata for Ubuntu 20.04
ARG TIMEZONE=Europe/Berlin
RUN export TZ=$TIMEZONE && echo $TZ > /etc/timezone && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    apt-get -yy update && apt-get -yy install wget tzdata && rm -rf /var/lib/apt/lists/*
# The version to fetch the package for: vX.Y.Z
ARG version
RUN echo "$version" | grep "v[0-9]\+\.[0-9]\+\.[0-9]\+" > /dev/null # The version must follow the format vX.Y.Z
RUN wget -q -O libprecice.deb https://github.com/precice/precice/releases/download/${version}/libprecice`echo ${version} | sed 's/v\([0-9]\+\)\.\([0-9]\+\.[0-9]\+\)/\1_\1.\2/'`_focal.deb && \
    ( dpkg -i libprecice.deb || apt-get install -fyy ) && \
    rm libprecice.deb && \
    binprecice xml > /dev/null # Make sure the installation is functional
