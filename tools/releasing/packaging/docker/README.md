# preCICE Dockerfiles

The point of these dockerfiles is to build a Ubuntu-based distribution of preCICE.

## Release images `release.dockerfile`

This Dockerfile uses named releases of preCICE such as `2.1.1` and installs the attached debian package in the container.

Use the following to build the `2.1.1` image locally:
```
cd tools/releasing/packaging/docker/
docker build -f release.dockerfile  -t precice/precice:2.1.1 --build-arg=version=2.1.1 .
```

## Nightly release images `nightly.dockerfile`

This Dockerfile builds the current develop version of preCICE and installs the debian package in the container.

```
cd tools/releasing/packaging/docker/
docker build -f nightly.dockerfile  -t precice/precice:nightly .
```
