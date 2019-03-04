#!/bin/sh
# Installs dependencies of preCICE
# boost, Eigen and PETSc
# if not already cached by Travis.

set -x
set -e

LOCAL_INSTALL=$1

# Don't test for $LOCAL_INSTALL, because it's created by the cacher.
if [ ! -d $LOCAL_INSTALL/include ]; then
    mkdir -p $LOCAL_INSTALL/include $LOCAL_INSTALL/lib

    # Download and extract Eigen
    mkdir -p $LOCAL_INSTALL/eigen3
    wget -nv http://bitbucket.org/eigen/eigen/get/3.3.2.tar.bz2 -O - | tar xj -C $LOCAL_INSTALL/eigen3 --strip-components=1 eigen-eigen-da9b4e14c255

    # Download, compile and install Boost
    wget -nv 'https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2' -O - | tar xj
    cd boost_1_65_1
    ./bootstrap.sh --prefix=$LOCAL_INSTALL > ~/boost.bootstrap
    ./b2 -j2 --with-program_options --with-test --with-filesystem --with-log install > ~/boost.b2

    # Download and compile PETSc
    cd $LOCAL_INSTALL
    git clone -b maint https://bitbucket.org/petsc/petsc petsc
    cd petsc
    export PETSC_ARCH=arch-linux2-c-debug
    python2 configure --with-debugging=1 --with-64-bit-indices > ~/petsc.configure
    make > ~/petsc.make
fi

# get modern cmake
if [ ! -d $LOCAL_INSTALL/cmake ]; then
    mkdir -p ${LOCAL_INSTALL}/cmake
    CMAKE_URL="http://www.cmake.org/files/v3.10/cmake-3.10.1-Linux-x86_64.tar.gz"
    mkdir -p ${LOCAL_INSTALL}/cmake
    wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C ${LOCAL_INSTALL}/cmake
fi
