# preCICE Change Log

All notable changes to this project will be documented in this file.

## develop

- Remove the `staticlib` and `bin` from the default targets to reduce the building time and storage requirements.
- Change the `-b` option to `-t` in the script `compileAndTest.py`.
- Change build types to mixed case, i.e. ```Debug``` and ```Release```. Old versions are retained for backward compatibility.
- Add experimental CMake build control files.

## 1.0.3
- Fix compilation for boost 1.66, see issue #93.

## 1.0.2
- Fix bug in the mesh repartitioning for plane-like coupling interfaces and small gaps between both sides.

## 1.0.1
- Fix compilation issue with the python interface.
