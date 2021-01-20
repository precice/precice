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
CACHE_CCACHE_TOKEN=$LOCAL_INSTALL/done-ccache

# Download and extract Eigen
if [ ! -f $CACHE_EIGEN_TOKEN ]; then
    # Cleanup
    rm -rf $LOCAL_INSTALL/eigen3
    mkdir $LOCAL_INSTALL/eigen3
    # Download
    wget -nv https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2 -O - | tar xj -C $LOCAL_INSTALL/eigen3 --strip-components=1
    # Create token
    touch $CACHE_EIGEN_TOKEN
fi

# Download, compile and install Boost
if [ ! -f $CACHE_BOOST_TOKEN ]; then
    # Cleanup
    rm -rf $LOCAL_INSTALL/boost
    mkdir $LOCAL_INSTALL/boost
    # Download
    wget -nv 'https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2' -O - | tar xj -C $LOCAL_INSTALL/boost
    # Compile and install
    cd $LOCAL_INSTALL/boost/boost_1_65_1
    ./bootstrap.sh --prefix=$LOCAL_INSTALL > ~/boost.bootstrap
    ./b2 -j2 --with-program_options --with-test --with-filesystem --with-log install > ~/boost.b2
    cd $LOCAL_INSTALL
    # Cleanup
    rm -rf $LOCAL_INSTALL/boost
    # Create token
    touch $CACHE_BOOST_TOKEN
fi

# Download and compile PETSc
if [ ! -f $CACHE_PETSC_TOKEN ]; then
    # Cleanup
    rm -rf $LOCAL_INSTALL/petsc
    mkdir $LOCAL_INSTALL/petsc
    # Download
    wget -nv 'https://gitlab.com/petsc/petsc/-/archive/maint/petsc-maint.tar.bz2' -O - | tar xj -C $LOCAL_INSTALL/petsc --strip-components=1
    # Configure and compile
    cd $LOCAL_INSTALL/petsc
    export PETSC_ARCH=arch-linux2-c-debug
    python configure --with-debugging=1 --with-64-bit-indices > ~/petsc.configure
    make > ~/petsc.make
    cd $LOCAL_INSTALL
    # Create token
    touch $CACHE_PETSC_TOKEN
fi

# Download CMake 3.10.2
if [ ! -f $CACHE_CMAKE_TOKEN ]; then
    # Cleanup
    rm -rf $LOCAL_INSTALL/cmake
    mkdir ${LOCAL_INSTALL}/cmake
    # Download
    CMAKE_URL="http://www.cmake.org/files/v3.10/cmake-3.10.2-Linux-x86_64.tar.gz"
    wget --no-check-certificate --quiet -O - ${CMAKE_URL} | tar --strip-components=1 -xz -C ${LOCAL_INSTALL}/cmake
    # Check version
    $LOCAL_INSTALL/cmake/bin/cmake --version
    # Create token
    touch $CACHE_CMAKE_TOKEN
fi


# Download and compile ccache 3.6
if [ ! -f $CACHE_CCACHE_TOKEN ]; then
    # Cleanup
    rm -rf $LOCAL_INSTALL/ccache
    mkdir ${LOCAL_INSTALL}/ccache
    # Download
    CCACHE_URL="https://www.samba.org/ftp/ccache/ccache-3.6.tar.gz"
    wget --no-check-certificate --quiet -O - ${CCACHE_URL} | tar --strip-components=1 -xz -C ${LOCAL_INSTALL}/ccache
    # Configure and compile
    cd $LOCAL_INSTALL/ccache
    ./configure --prefix=$LOCAL_INSTALL --disable-man
    make -j $(nproc)
    make install
    strip $LOCAL_INSTALL/bin/ccache
    cd $LOCAL_INSTALL
    # Cleanup
    rm -rf $LOCAL_INSTALL/ccache
    # Check version
    $LOCAL_INSTALL/bin/ccache --version
    # Create token
    touch $CACHE_CCACHE_TOKEN
fi
