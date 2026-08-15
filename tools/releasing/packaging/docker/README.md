# preCICE Dockerfiles

The point of these dockerfiles is to build a Ubuntu-based distribution of preCICE.

All images contain the user `precice`, allowing them to run executables using MPI.

## Release images `release.dockerfile`

This Dockerfile uses the matching debian package of preCICE from the build context and installs the debian package in the container.

Use the following to build the `3.1.2` image locally:
```
gh release download v3.1.2
docker build -f release.dockerfile  -t precice/precice:3.1.2 .
```

To upgrade the Ubuntu version, change the version of the baseimage `FROM ubuntu` and the `CODENAME` (output of `lsb_release -sc)`.

## Nightly release images `nightly.dockerfile`

This Dockerfile builds the current develop version of preCICE and installs the debian package in the container.

```
cd tools/releasing/packaging/docker/
docker build -f nightly.dockerfile  -t precice/precice:nightly .
```
