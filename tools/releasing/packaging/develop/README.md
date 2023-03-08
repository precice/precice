# preCICE develop Dockerfile

The point of this dockerfile is to build a Ubuntu-based test distribution of the current preCICE develop.

Use the following to build the image locally:
```
cd tools/releasing/packaging/develop/
docker build -t precice/precice:develop .
```
