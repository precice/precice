# preCICE Dockerfile

The point of this dockerfile is to build a Ubuntu-based distribution of preCICE.
It uses named releases of precice such as `v2.1.1` and installs the attached debian package in the container.

Use the following to build the `v2.1.1` image locally:
```
cd tools/releasing/packaging/docker/
docker build -t precice/precice --build-arg=version=v2.1.1 .
```
