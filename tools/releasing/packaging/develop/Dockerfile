# Dockerfile to build a ubuntu image containing the installed develop

FROM precice/ci-ubuntu-2204
RUN git clone --depth=1 --branch=develop https://github.com/precice/precice.git precice-build && \
    cd precice-build && \
    cmake -DBUILD_TESTING=OFF . && \
    make -j $(nproc) install && \
    cd .. && \
    rm -r precice-build
