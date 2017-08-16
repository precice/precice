#!/bin/sh
# Installs dependencies of preCICE
# boost and Eigen
# if not already cached by Travis.

set -x
set -e

LOCAL_INSTALL=$1

# Don't test for $LOCAL_INSTALL, because it's created by the cacher.
if [ ! -d $LOCAL_INSTALL/include ]; then
    mkdir -p $LOCAL_INSTALL/include $LOCAL_INSTALL/lib

    # Download and extract Eigen
    wget http://bitbucket.org/eigen/eigen/get/3.3.2.tar.bz2 -O - | tar xj -C $LOCAL_INSTALL/include --strip-components=1 eigen-eigen-da9b4e14c255/Eigen

    # Download, compile and install Boost
    wget 'http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.bz2' -O - | tar xj
    cd boost_1_60_0
    ./bootstrap.sh --prefix=$LOCAL_INSTALL > ~/boost.bootstrap
    ./b2 -j2 --with-program_options --with-test --with-filesystem --with-log install > ~/boost.b2
fi

