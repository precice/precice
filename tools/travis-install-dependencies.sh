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
    wget -nv http://bitbucket.org/eigen/eigen/get/3.3.2.tar.bz2 -O - | tar xj -C $LOCAL_INSTALL/include --strip-components=1 eigen-eigen-da9b4e14c255/Eigen

    # Download, compile and install Boost
    wget -nv 'https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2' -O - | tar xj
    cd boost_1_65_1
    ./bootstrap.sh --prefix=$LOCAL_INSTALL > ~/boost.bootstrap
    ./b2 -j2 --with-program_options --with-test --with-filesystem --with-log --with-iostreams install > ~/boost.b2

    # Download and compile PETSc
    cd $LOCAL_INSTALL
    git clone -b maint https://bitbucket.org/petsc/petsc petsc
    cd petsc
    export PETSC_ARCH=arch-linux2-c-debug
    python2 configure --with-debugging=1 --with-64-bit-indices > ~/petsc.configure
    make > ~/petsc.make
fi

