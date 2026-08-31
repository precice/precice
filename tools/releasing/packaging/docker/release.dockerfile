# Dockerfile to build a ubuntu image containing the installed Debian package of a release
#
# The build context must include the debian packages of the release
# Use gh release download vX.Y.Z

# Update the following two lines when bumping the base version.
FROM ubuntu:24.04
ENV CODENAME=noble

# Add the precice user
RUN useradd -m -s /bin/bash precice

# Copy all debian packages from the build context.
COPY ./libprecice*.deb /debs/

# Fix the installation of tzdata for Ubuntu
ARG TIMEZONE=Europe/Berlin
ENV TZ=${TIMEZONE}

# Install tzdata and the matching debian package and cleanup
RUN echo ${TZ} > /etc/timezone && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    apt-get -yy update && apt-get -yy install tzdata wget && \
    apt-get -yy install /debs/libprecice*_${CODENAME}.deb && \
    rm -rf /var/lib/apt/lists/* /debs

# Make sure the installation is functional
RUN precice-version
