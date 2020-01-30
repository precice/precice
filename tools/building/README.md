# Tools for assisting with building preCICE

This directory contains optional tools that are useful for developers
building preCICE often.

- `compileAndTest.py`: Build and test preCICE automatically.
  Useful, for example, for identifying the source of a regression (using `git blame`)
- `updateSourceFiles.py`: Update the list of source files that CMake needs to build.
  Useful when adding/removing source files.