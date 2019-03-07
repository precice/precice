#!/bin/sh
# Installs dependencies of preCICE
# boost, Eigen and PETSc
# if not already cached by Travis.

set -x
set -e

LOCAL_INSTALL=$1

CACHE_EIGEN_TOKEN=$LOCAL_INSTALL/done-eigen3
CACHE_BOOST_TOKEN=$LOCAL_INSTALL/done-boost
CACHE_PETSC_TOKEN=$LOCAL_INSTALL/done-petsc
CACHE_CMAKE_TOKEN=$LOCAL_INSTALL/done-cmake

if [ -f $CACHE_EIGEN_TOKEN ]; then
    echo "[HIT] eigen3"
else
    echo "[MISS] eigen3"
    # Download and extract Eigen
    rm -r $LOCAL_INSTALL/eigen3
    mkdir $LOCAL_INSTALL/eigen3
    wget -nv http://bitbucket.org/eigen/eigen/get/3.3.2.tar.bz2 -O - | tar xj -C $LOCAL_INSTALL/eigen3 --strip-components=1 eigen-eigen-da9b4e14c255
    ls -hl $LOCAL_INSTALL/eigen3
    touch $CACHE_EIGEN_TOKEN
fi

# Download, compile and install Boost
if [ -f $CACHE_BOOST_TOKEN ]; then
    echo "[HIT] boost"
else
    echo "[MISS] boost"
    rm -r $LOCAL_INSTALL/boost
    wget -nv 'https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2' -O - | tar xj -C $LOCAL_INSTALL/boost
    cd $LOCAL_INSTALL/boost/boost_1_65_1
    ./bootstrap.sh --prefix=$LOCAL_INSTALL > ~/boost.bootstrap
    ./b2 -j2 --with-program_options --with-test --with-filesystem --with-log install > ~/boost.b2
    # rm -r $LOCAL_INSTALL/boost
    touch $CACHE_BOOST_TOKEN
fi

# Download and compile PETSc
if [ -f $CACHE_PETSC_TOKEN ]; then
    echo "[HIT] petsc"
else
    echo "[MISS] petsc"
    rm -r $LOCAL_INSTALL/petsc
    git clone -b maint https://bitbucket.org/petsc/petsc $LOCAL_INSTALL/petsc
    cd $LOCAL_INSTALL/petsc
    export PETSC_ARCH=arch-linux2-c-debug
    python2 configure --with-debugging=1 --with-64-bit-indices > ~/petsc.configure
    make > ~/petsc.make
    # rm -r $LOCAL_INSTALL/petsc
    touch $CACHE_PETSC_TOKEN
fi

# Download CMake 3.10.1
if [ -f $CACHE_CMAKE_TOKEN ]; then
    echo "[HIT] cmake"
else
    echo "[MISS] cmake"
    rm -r $LOCAL_INSTALL/cmake
    CMAKE_URL="http://www.cmake.org/files/v3.10/cmake-3.10.1-Linux-x86_64.tar.gz"
    mkdir ${LOCAL_INSTALL}/cmake
    wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C ${LOCAL_INSTALL}/cmake
    $LOCAL_INSTALL/cmake/bin/cmake --version
    touch $CACHE_CMAKE_TOKEN
fi
